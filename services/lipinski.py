"""
app/services/lipinski.py
========================
Lipinski's Rule of Five (Ro5) evaluator.

Background
----------
Proposed by Christopher Lipinski at Pfizer in 1997, the Rule of Five
describes physicochemical properties associated with good oral
bioavailability in humans.  The "five" in the name refers to the fact
that all threshold values are multiples of five:

  Rule 1 — Molecular Weight      ≤ 500 Da
  Rule 2 — LogP (Wildman-Crippen) ≤ 5
  Rule 3 — H-bond Donors         ≤ 5
  Rule 4 — H-bond Acceptors      ≤ 10

A molecule with **≤ 1 violation** is generally considered orally drug-like.
Molecules with **≥ 2 violations** face significant ADME hurdles for oral dosing,
though natural products and macrolides are common exceptions.

Note: The Ro5 predicts *passive intestinal absorption*, not bioavailability
per se.  Active transport, first-pass metabolism, and formulation factors
are not captured.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class LipinskiResult:
    """
    Detailed outcome of a Lipinski Ro5 evaluation.

    Each rule is recorded individually to support granular reporting
    and downstream lead-optimisation guidance.
    """

    violations: int
    rule_results: Dict[str, bool]    # True = rule PASSED (within limit)
    violation_details: List[str]     # Human-readable explanation of each failure


# Thresholds are module-level constants so they can be referenced in tests
MW_LIMIT: float = 500.0
LOGP_LIMIT: float = 5.0
HBD_LIMIT: int = 5
HBA_LIMIT: int = 10


def evaluate_lipinski(
    molecular_weight: float,
    logp: float,
    hbond_donors: int,
    hbond_acceptors: int,
) -> LipinskiResult:
    """
    Apply Lipinski's Rule of Five to a set of pre-computed descriptors.

    Parameters
    ----------
    molecular_weight:
        Exact molecular weight (g/mol).
    logp:
        Wildman-Crippen log P.
    hbond_donors:
        Number of NH + OH groups.
    hbond_acceptors:
        Number of N + O atoms (Ertl definition).

    Returns
    -------
    LipinskiResult
        Structured outcome with per-rule pass/fail flags and violation narratives.
    """
    rule_results: Dict[str, bool] = {}
    violation_details: List[str] = []

    # ------------------------------------------------------------------
    # Rule 1 — Molecular Weight
    # ------------------------------------------------------------------
    mw_pass = molecular_weight <= MW_LIMIT
    rule_results["MW ≤ 500 Da"] = mw_pass
    if not mw_pass:
        violation_details.append(
            f"MW = {molecular_weight:.2f} Da exceeds the 500 Da threshold. "
            "High molecular weight reduces passive membrane permeability and "
            "increases the likelihood of P-gp efflux."
        )

    # ------------------------------------------------------------------
    # Rule 2 — LogP
    # ------------------------------------------------------------------
    logp_pass = logp <= LOGP_LIMIT
    rule_results["LogP ≤ 5"] = logp_pass
    if not logp_pass:
        violation_details.append(
            f"LogP = {logp:.2f} exceeds the upper limit of 5. "
            "High lipophilicity is associated with poor aqueous solubility, "
            "elevated non-specific plasma protein binding, and increased "
            "risk of CYP inhibition and hepatotoxicity."
        )

    # ------------------------------------------------------------------
    # Rule 3 — H-bond Donors
    # ------------------------------------------------------------------
    hbd_pass = hbond_donors <= HBD_LIMIT
    rule_results["HBD ≤ 5"] = hbd_pass
    if not hbd_pass:
        violation_details.append(
            f"H-bond donors = {hbond_donors} exceeds the limit of 5. "
            "Excess hydrogen-bond donors penalise membrane desolvation energy, "
            "reducing passive transcellular permeability."
        )

    # ------------------------------------------------------------------
    # Rule 4 — H-bond Acceptors
    # ------------------------------------------------------------------
    hba_pass = hbond_acceptors <= HBA_LIMIT
    rule_results["HBA ≤ 10"] = hba_pass
    if not hba_pass:
        violation_details.append(
            f"H-bond acceptors = {hbond_acceptors} exceeds the limit of 10. "
            "High acceptor count increases hydrophilicity beyond the optimal "
            "range for passive absorption, favouring aqueous solvation over "
            "membrane partitioning."
        )

    violations = sum(1 for passed in rule_results.values() if not passed)

    return LipinskiResult(
        violations=violations,
        rule_results=rule_results,
        violation_details=violation_details,
    )


def is_drug_like(result: LipinskiResult) -> bool:
    """Return True if the molecule satisfies the Ro5 drug-likeness criterion (≤ 1 violation)."""
    return result.violations <= 1
