"""
app/services/medicinal_reasoning.py
=====================================
Medicinal chemistry reasoning engine.

This is the intelligence layer of MolecularMind.  It transforms raw
numerical descriptors into structured, expert-level pharmacokinetic
interpretations.

Design philosophy
-----------------
Real medicinal chemists think in terms of competing optimisation pressures:
improving LogP to boost membrane permeability invariably worsens aqueous
solubility; reducing MW to aid absorption may remove pharmacophore atoms
needed for target binding.  This engine makes those tensions explicit.

The reasoning functions follow a threshold-based decision tree derived from
landmark publications:

  * Lipinski et al., Adv Drug Deliv Rev (1997, 2001) — Ro5
  * Veber et al., J Med Chem (2002) — rotatable bonds / TPSA
  * Egan et al., J Med Chem (2000) — absorption / permeability ellipse
  * Kelder et al., Pharm Res (1999) — BBB penetration TPSA < 60–90 Å²
  * Leeson & Springthorpe, Nat Rev Drug Discov (2007) — lipophilicity efficiency
  * Lovering et al., J Med Chem (2009) — Fsp³ and saturation

Each function is independent and returns a plain string, making it
straightforward to unit-test in isolation.
"""

from __future__ import annotations

import logging
from typing import List

from app.services.descriptor_calculator import RawDescriptors
from app.services.lipinski import LipinskiResult, is_drug_like

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_analysis(descriptors: RawDescriptors, lipinski: LipinskiResult) -> dict:
    """
    Orchestrate the full medicinal chemistry reasoning pipeline.

    Parameters
    ----------
    descriptors:
        Computed physicochemical descriptor profile.
    lipinski:
        Result of the Lipinski Ro5 evaluation.

    Returns
    -------
    dict
        Keys match the MedicinalChemistryAnalysis Pydantic schema.
    """
    conflicts = _identify_optimisation_conflicts(descriptors, lipinski)

    return {
        "absorption": _reason_absorption(descriptors),
        "permeability": _reason_permeability(descriptors),
        "solubility": _reason_solubility(descriptors),
        "bbb_penetration": _reason_bbb(descriptors),
        "drug_likeness": _reason_drug_likeness(descriptors, lipinski),
        "flexibility_risk": _reason_flexibility(descriptors),
        "optimization_conflicts": conflicts,
        "overall_verdict": _overall_verdict(descriptors, lipinski, conflicts),
    }


# ---------------------------------------------------------------------------
# Individual reasoning functions
# ---------------------------------------------------------------------------

def _reason_absorption(d: RawDescriptors) -> str:
    """
    Predict oral absorption using MW, TPSA, and rotatable bonds.

    Reference: Veber et al. (2002) — TPSA ≤ 140 Å² and rotatable bonds ≤ 10
    are better predictors of oral bioavailability in rats than MW alone.
    """
    flags: List[str] = []

    # MW contribution
    if d.molecular_weight <= 300:
        mw_comment = (
            f"The low molecular weight ({d.molecular_weight:.2f} Da) provides an "
            "excellent foundation for passive intestinal absorption."
        )
    elif d.molecular_weight <= 500:
        mw_comment = (
            f"Molecular weight of {d.molecular_weight:.2f} Da is within the Ro5 "
            "acceptable range for oral absorption."
        )
    elif d.molecular_weight <= 700:
        mw_comment = (
            f"Elevated molecular weight ({d.molecular_weight:.2f} Da) exceeds the "
            "Lipinski threshold (500 Da), suggesting reduced passive permeability. "
            "Dose escalation or prodrug strategies may be warranted."
        )
    else:
        mw_comment = (
            f"Very high molecular weight ({d.molecular_weight:.2f} Da) is strongly "
            "unfavourable for oral bioavailability. Reformulation as a biologic or "
            "non-oral route of administration should be considered."
        )
    flags.append(mw_comment)

    # TPSA contribution (Veber rule)
    if d.tpsa <= 60:
        flags.append(
            f"TPSA of {d.tpsa:.1f} Å² is low — excellent for intestinal epithelial "
            "permeation via passive transcellular transport."
        )
    elif d.tpsa <= 90:
        flags.append(
            f"TPSA of {d.tpsa:.1f} Å² is in the optimal oral-absorption window (< 90 Å²)."
        )
    elif d.tpsa <= 140:
        flags.append(
            f"TPSA of {d.tpsa:.1f} Å² approaches the Veber limit (140 Å²). "
            "Intestinal absorption may be reduced by 30–50% compared to low-TPSA analogues."
        )
    else:
        flags.append(
            f"TPSA of {d.tpsa:.1f} Å² exceeds 140 Å² — passive oral absorption is "
            "predicted to be very poor. The molecule is likely to remain in the GI lumen "
            "or require an active transporter for uptake."
        )

    # Rotatable bonds
    if d.rotatable_bonds > 10:
        flags.append(
            f"High rotatable bond count ({d.rotatable_bonds}) exceeds Veber's threshold "
            "of 10, increasing the conformational entropy penalty upon membrane desolvation."
        )
    elif d.rotatable_bonds > 7:
        flags.append(
            f"Moderate rotatable bond count ({d.rotatable_bonds}) — approaching Veber's limit. "
            "Monitor carefully during analogue synthesis."
        )

    return " ".join(flags)


def _reason_permeability(d: RawDescriptors) -> str:
    """
    Assess membrane permeability using the Egan 'egg' model parameters:
    LogP (x-axis) and TPSA (y-axis).

    Optimal Egan space: −1 ≤ LogP ≤ 5.88, TPSA ≤ 131.6 Å²
    """
    parts: List[str] = []

    if d.logp < 0:
        parts.append(
            f"LogP = {d.logp:.2f} is strongly negative, indicating high hydrophilicity. "
            "The molecule is poorly suited for passive lipid-bilayer permeation. "
            "Uptake will depend on paracellular transport or active influx transporters."
        )
    elif d.logp <= 2:
        parts.append(
            f"LogP = {d.logp:.2f} is in the hydrophilic range. Passive transcellular "
            "permeability is moderate; paracellular routes may contribute."
        )
    elif d.logp <= 3.5:
        parts.append(
            f"LogP = {d.logp:.2f} sits in the ideal window for membrane permeability. "
            "Adequate lipophilicity for passive transcellular diffusion without "
            "excessive non-specific binding or solubility penalties."
        )
    elif d.logp <= 5:
        parts.append(
            f"LogP = {d.logp:.2f} is moderately lipophilic — favours membrane "
            "permeability but increases risk of non-specific protein binding and "
            "metabolic liability (CYP3A4 in particular)."
        )
    else:
        parts.append(
            f"LogP = {d.logp:.2f} is high — although lipid membrane partitioning "
            "is favoured, aqueous solubility will be severely limited, creating a "
            "dissolution bottleneck that offsets the permeability advantage. "
            "Risk of non-specific toxicity and hERG channel binding increases."
        )

    # TPSA contribution to permeability
    if d.tpsa > 120:
        parts.append(
            "The high TPSA compounds the permeability challenge: polar surface area "
            "in excess of 120 Å² strongly disfavours passive transcellular transport."
        )

    return " ".join(parts)


def _reason_solubility(d: RawDescriptors) -> str:
    """
    Estimate aqueous solubility risk from LogP and TPSA.

    General heuristics:
      LogP < 2    → good solubility (though ionisation state matters)
      LogP 2–4    → moderate solubility
      LogP > 4    → poor solubility risk
      TPSA > 80   → increased solvation, improved solubility

    Note: ESOL (Delaney 2004) is the gold-standard computational model;
    these rules are first-order approximations.
    """
    parts: List[str] = []

    if d.logp < 1:
        parts.append(
            f"LogP = {d.logp:.2f} predicts good intrinsic aqueous solubility. "
            "Solubility is unlikely to limit bioavailability for this compound."
        )
    elif d.logp <= 2:
        parts.append(
            f"LogP = {d.logp:.2f} suggests adequate solubility for most oral formulations."
        )
    elif d.logp <= 4:
        parts.append(
            f"LogP = {d.logp:.2f} indicates moderate solubility. Aqueous solubility should "
            "be confirmed experimentally (miniaturised nephelometry or thermodynamic solubility assay). "
            "Formulation optimisation (amorphous solid dispersion, co-solvent, surfactant) "
            "may be needed for high-dose scenarios."
        )
    else:
        parts.append(
            f"LogP = {d.logp:.2f} predicts poor aqueous solubility — a major liability "
            "for oral dosing. High lipophilicity drives precipitation in the GI aqueous "
            "environment. Biopharmaceutical Classification System (BCS) Class II or IV "
            "behaviour is probable."
        )

    # TPSA as a counter-balancing solubility driver
    if d.tpsa >= 80:
        parts.append(
            f"A TPSA of {d.tpsa:.1f} Å² provides compensating hydrophilic surface area "
            "that partially mitigates the lipophilicity-driven solubility concern."
        )
    elif d.tpsa < 40 and d.logp > 3:
        parts.append(
            "Very low TPSA combined with elevated LogP creates a particularly unfavourable "
            "solubility profile — prioritise structural modifications that introduce "
            "polar, solubilising groups (e.g., secondary amines, hydroxyl groups)."
        )

    return " ".join(parts)


def _reason_bbb(d: RawDescriptors) -> str:
    """
    Estimate blood-brain barrier (BBB) penetration likelihood.

    CNS penetration heuristics (Kelder 1999, Pajouhesh & Lenz 2005):
      TPSA < 60 Å²  → high BBB penetration probability
      TPSA < 90 Å²  → moderate probability
      TPSA > 90 Å²  → low probability
      MW < 400      → favoured for CNS
      −0.5 < LogP < 5 → CNS window
      HBD ≤ 3       → supports CNS penetration
    """
    parts: List[str] = []

    # TPSA — primary determinant
    if d.tpsa < 60:
        parts.append(
            f"TPSA = {d.tpsa:.1f} Å² (< 60 Å² threshold) strongly supports CNS penetration "
            "via passive transcellular diffusion across the blood-brain barrier."
        )
    elif d.tpsa < 90:
        parts.append(
            f"TPSA = {d.tpsa:.1f} Å² allows moderate BBB penetration (60–90 Å² borderline zone). "
            "CNS bioavailability may be sufficient for central targets with adequate potency."
        )
    else:
        parts.append(
            f"TPSA = {d.tpsa:.1f} Å² (> 90 Å²) strongly disfavours BBB penetration. "
            "If CNS activity is required, polar group removal or bioisosteric replacement "
            "strategies (e.g., replacing COOH with tetrazole) should be explored."
        )

    # MW for CNS
    if d.molecular_weight <= 400:
        parts.append(
            f"MW = {d.molecular_weight:.2f} Da is within the preferred CNS molecular weight range (< 400 Da)."
        )
    elif d.molecular_weight <= 500:
        parts.append(
            f"MW = {d.molecular_weight:.2f} Da is marginally elevated for CNS penetration — "
            "the 400 Da soft limit reflects tighter efflux constraints at the BBB vs gut epithelium."
        )
    else:
        parts.append(
            f"MW = {d.molecular_weight:.2f} Da exceeds both the Ro5 and CNS MW limits, "
            "substantially diminishing BBB penetration probability."
        )

    # LogP CNS window
    if 1 <= d.logp <= 3:
        parts.append("LogP is within the optimal CNS lipophilicity window (1–3).")
    elif d.logp > 5:
        parts.append(
            "High LogP may favour BBB membrane partitioning in theory but will "
            "increase P-gp efflux liability and non-specific CNS tissue binding, "
            "reducing free brain concentration."
        )

    return " ".join(parts)


def _reason_drug_likeness(d: RawDescriptors, lipinski: LipinskiResult) -> str:
    """
    Overall drug-likeness assessment integrating Ro5, Veber, Ghose, and Fsp³ criteria.
    """
    parts: List[str] = []

    # Lipinski Ro5
    if lipinski.violations == 0:
        parts.append(
            "The molecule satisfies all four Lipinski Ro5 criteria with zero violations — "
            "an excellent starting point for an orally bioavailable small molecule drug."
        )
    elif lipinski.violations == 1:
        parts.append(
            f"One Lipinski Ro5 violation detected. Per the original rule, one violation "
            "is permissible (particularly for substrates of active transporters). "
            f"Violated parameter(s): {[k for k, v in lipinski.rule_results.items() if not v]}."
        )
    else:
        part = (
            f"{lipinski.violations} Lipinski Ro5 violations detected. "
            "Oral bioavailability via passive absorption is predicted to be poor. "
            f"Violated rule(s): {[k for k, v in lipinski.rule_results.items() if not v]}. "
        )
        if lipinski.violations >= 3:
            part += (
                "Consider non-oral delivery routes, prodrug approaches, "
                "or fundamental scaffold redesign."
            )
        parts.append(part)

    # Ghose filter (MW 160–480, LogP −0.4–5.6, MR 40–130, atom count 20–70)
    ghose_issues: List[str] = []
    if not (160 <= d.molecular_weight <= 480):
        ghose_issues.append(f"MW ({d.molecular_weight:.1f} Da) outside 160–480 Da")
    if not (-0.4 <= d.logp <= 5.6):
        ghose_issues.append(f"LogP ({d.logp:.2f}) outside −0.4–5.6")
    if not (40 <= d.molar_refractivity <= 130):
        ghose_issues.append(f"MR ({d.molar_refractivity:.1f}) outside 40–130")
    if ghose_issues:
        parts.append(
            "Ghose filter deviations noted: " + "; ".join(ghose_issues) + ". "
            "The Ghose filter is more stringent than Ro5 — deviations here do not "
            "disqualify a candidate but warrant monitoring."
        )

    # Fsp³ — molecular saturation
    if d.fraction_csp3 >= 0.25:
        parts.append(
            f"Fsp³ = {d.fraction_csp3:.2f} (≥ 0.25) indicates good molecular saturation. "
            "Higher Fsp³ correlates with reduced flat-aromatic character, "
            "improved selectivity, and lower CYP inhibition risk (Lovering et al., 2009)."
        )
    elif d.fraction_csp3 >= 0.1:
        parts.append(
            f"Fsp³ = {d.fraction_csp3:.2f} is in the intermediate range. "
            "The molecule has moderate aromatic character — consider sp³-enriching "
            "modifications to improve selectivity and developability."
        )
    else:
        parts.append(
            f"Low Fsp³ = {d.fraction_csp3:.2f} indicates a highly flat/aromatic scaffold. "
            "Such compounds are associated with elevated CYP inhibition, poor crystal packing, "
            "lower clinical success rates, and reduced aqueous solubility."
        )

    return " ".join(parts)


def _reason_flexibility(d: RawDescriptors) -> str:
    """
    Assess conformational flexibility risk from rotatable bond count.

    Veber et al. (2002): ≤ 10 rotatable bonds for good oral bioavailability.
    Congreve et al. (2003): "Rule of Three" for fragments uses ≤ 3.
    High rotatable bonds → high conformational entropy penalty upon binding →
    reduced binding affinity / poor enthalpy-entropy compensation.
    """
    if d.rotatable_bonds <= 3:
        return (
            f"Rotatable bonds = {d.rotatable_bonds} — very rigid scaffold. "
            "Low conformational entropy penalty upon protein binding. "
            "Pre-organisation can translate to enthalpically favourable binding, "
            "but over-rigidity may prevent induced-fit and reduce binding breadth."
        )
    elif d.rotatable_bonds <= 7:
        return (
            f"Rotatable bonds = {d.rotatable_bonds} — optimal flexibility for "
            "a drug-like molecule. Sufficient conformational freedom for "
            "binding-pocket adaptation without excessive entropy cost."
        )
    elif d.rotatable_bonds <= 10:
        return (
            f"Rotatable bonds = {d.rotatable_bonds} — approaching the Veber limit (10). "
            "Oral bioavailability is unlikely to be compromised, but binding entropy "
            "penalty is becoming significant. Monitor carefully."
        )
    elif d.rotatable_bonds <= 15:
        return (
            f"Rotatable bonds = {d.rotatable_bonds} exceeds Veber's guideline. "
            "Conformational flexibility is high — the molecule will pay a substantial "
            "entropy penalty upon adopting a binding-competent conformation. "
            "Macrocyclisation or intramolecular H-bond formation may rigidify the scaffold."
        )
    else:
        return (
            f"High rotatable bond count ({d.rotatable_bonds}) is a significant liability. "
            "The molecule resembles a 'flexible linker' with limited pre-organisation. "
            "Binding affinity will be entropically penalised, and oral absorption is "
            "predicted to be poor (Veber et al., 2002)."
        )


def _identify_optimisation_conflicts(
    d: RawDescriptors, lipinski: LipinskiResult
) -> List[str]:
    """
    Identify lead-optimisation tensions — competing pressures where improving
    one property necessarily worsens another.

    These are the most actionable outputs for a medicinal chemist.
    """
    conflicts: List[str] = []

    # Classic solubility vs. permeability trade-off
    if d.logp > 3 and d.tpsa < 60:
        conflicts.append(
            "Permeability-solubility conflict: High LogP favours membrane permeability "
            "but drives down aqueous solubility. Introducing a polar, ionisable group "
            "(pKa-matched to GI pH) can improve solubility without collapsing permeability."
        )

    if d.logp < 1 and d.tpsa > 100:
        conflicts.append(
            "Absorption conflict: Both low LogP and high TPSA independently reduce "
            "passive membrane permeability. The compound may rely entirely on active "
            "transport. Risk of inter-individual variability due to transporter polymorphisms."
        )

    # BBB vs. periphery selectivity
    if d.tpsa < 60 and lipinski.violations >= 1:
        conflicts.append(
            "CNS penetration vs. drug-likeness conflict: TPSA favours BBB penetration, "
            "but Ro5 violations compromise oral bioavailability. Improving Ro5 compliance "
            "(e.g., MW reduction) may require removing polar atoms that also support target binding."
        )

    # Rigidity vs. binding breadth
    if d.rotatable_bonds <= 3 and d.ring_count >= 4:
        conflicts.append(
            "Rigidity vs. selectivity: Highly rigid polycyclic scaffold may achieve "
            "excellent binding affinity for one target but reduce conformational adaptability, "
            "limiting broader target engagement or off-target avoidance."
        )

    # Lipophilicity efficiency
    if d.logp > 4 and d.molecular_weight > 450:
        conflicts.append(
            "Lipophilicity efficiency (LipE / LLE) concern: High LogP AND high MW jointly "
            "erode lipophilic efficiency (pIC50 − LogP). Potency is expensive in molecular "
            "property space — prioritise potency-per-LogP-unit in SAR optimisation."
        )

    # Fsp3 vs. permeability
    if d.fraction_csp3 > 0.5 and d.tpsa > 90:
        conflicts.append(
            "Saturation vs. polarity conflict: High Fsp³ (good for selectivity) "
            "combined with high TPSA creates a compound with limited permeability "
            "despite good sp³ character. A constrained ring strategy may maintain "
            "Fsp³ while reducing the polar surface area."
        )

    # Aromatic rings — pan-assay interference (PAINS) proxy
    if d.aromatic_rings >= 4:
        conflicts.append(
            "High aromatic ring count (≥ 4) is a pan-assay interference (PAINS) "
            "risk factor and correlates with reduced clinical success rates. "
            "Consider replacing one aromatic ring with an aliphatic or partially "
            "saturated analogue to improve the Fsp³ profile."
        )

    return conflicts


def _overall_verdict(
    d: RawDescriptors, lipinski: LipinskiResult, conflicts: List[str]
) -> str:
    """
    Synthesise a concise executive-level summary for medicinal chemistry decision-making.
    """
    drug_like = is_drug_like(lipinski)
    conflict_count = len(conflicts)

    if drug_like and conflict_count == 0:
        quality = "excellent"
        action = (
            "This compound is a strong hit-to-lead candidate. Proceed to in vitro ADME "
            "profiling (Caco-2/MDCK permeability, microsomal stability, plasma protein binding) "
            "and selectivity assessment."
        )
    elif drug_like and conflict_count <= 2:
        quality = "good"
        action = (
            "Drug-like profile with manageable optimisation tensions. "
            f"{conflict_count} competing physicochemical pressure(s) identified — "
            "address via targeted analogue synthesis guided by the conflict analysis above."
        )
    elif drug_like and conflict_count > 2:
        quality = "moderate"
        action = (
            "Drug-like by Ro5 but with multiple competing optimisation pressures. "
            "Prioritise the 1–2 most critical conflicts identified above before "
            "expanding the analogue series."
        )
    elif not drug_like and lipinski.violations == 2:
        quality = "borderline"
        action = (
            "Two Ro5 violations place this molecule in a challenging developability space. "
            "Evaluate whether the violations are intrinsic to the pharmacophore or "
            "removable through bioisosteric replacement. Consider formulation strategies."
        )
    else:
        quality = "poor"
        action = (
            f"With {lipinski.violations} Ro5 violations and {conflict_count} optimisation "
            "conflicts, this scaffold faces significant ADME barriers. "
            "Fundamental redesign or a non-oral delivery route is recommended."
        )

    return (
        f"Overall drug-likeness: {quality.upper()}. "
        f"MW = {d.molecular_weight:.1f} Da | LogP = {d.logp:.2f} | "
        f"TPSA = {d.tpsa:.1f} Å² | HBD/HBA = {d.hbond_donors}/{d.hbond_acceptors} | "
        f"RotBonds = {d.rotatable_bonds} | Fsp³ = {d.fraction_csp3:.2f} | "
        f"Ro5 violations = {lipinski.violations}. {action}"
    )
