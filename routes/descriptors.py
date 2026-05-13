"""
app/routes/descriptors.py
=========================
FastAPI router — Molecular Descriptor Evaluation endpoints.

This module contains **only** routing and orchestration logic.
All domain work is delegated to the service layer.  The router's
responsibility is:
  1. Accept and validate the request payload.
  2. Invoke services in the correct order.
  3. Assemble and return the response model.
  4. Handle errors with meaningful HTTP semantics.
"""

from __future__ import annotations

import logging

from fastapi import APIRouter, HTTPException, status

from app.config import settings
from app.models.schemas import (
    MedicinalChemistryAnalysis,
    MolecularDescriptors,
    MoleculeEvaluationResponse,
    SMILESRequest,
)
from app.services.descriptor_calculator import compute_descriptors
from app.services.lipinski import evaluate_lipinski
from app.services.medicinal_reasoning import build_analysis
from app.services.rdkit_engine import parse_smiles
from app.utils.validators import check_molecule_complexity, validate_smiles_string

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# POST /evaluate
# ---------------------------------------------------------------------------

@router.post(
    "/evaluate",
    response_model=MoleculeEvaluationResponse,
    status_code=status.HTTP_200_OK,
    summary="Evaluate a molecule from a SMILES string",
    description=(
        "Submit a SMILES string to receive a complete physicochemical descriptor "
        "profile and medicinal chemistry interpretation. "
        "The response includes Lipinski Ro5 evaluation, ADME prediction reasoning, "
        "BBB penetration likelihood, and lead-optimisation conflict analysis."
    ),
    responses={
        200: {"description": "Successful evaluation (valid or invalid molecule)"},
        422: {"description": "Request payload validation error (malformed JSON, missing field)"},
        500: {"description": "Unexpected server-side error"},
    },
)
async def evaluate_molecule(payload: SMILESRequest) -> MoleculeEvaluationResponse:
    """
    Full molecular evaluation pipeline:

    1. Pre-validation (character set, bracket balance, InChI detection)
    2. RDKit parsing and sanitisation
    3. Molecule complexity checks (atom count, fragments)
    4. Descriptor computation
    5. Lipinski Ro5 evaluation
    6. Medicinal chemistry reasoning
    7. Response assembly
    """
    smiles = payload.smiles
    logger.info("Evaluation request received — SMILES: %s", smiles[:120])

    # -----------------------------------------------------------------------
    # Step 1 — Pre-parse string validation
    # -----------------------------------------------------------------------
    pre_check = validate_smiles_string(smiles, max_atom_hint=settings.MAX_ATOM_COUNT)

    if not pre_check.is_valid:
        logger.warning("Pre-validation failed: %s", pre_check.errors)
        return MoleculeEvaluationResponse(
            smiles=smiles,
            valid_molecule=False,
            error_detail=" | ".join(pre_check.errors),
            warnings=pre_check.warnings,
        )

    # Use the sanitised (whitespace-stripped) SMILES from now on
    clean_smiles = pre_check.sanitised_smiles or smiles

    # -----------------------------------------------------------------------
    # Step 2 — RDKit parsing
    # -----------------------------------------------------------------------
    parse_result = parse_smiles(clean_smiles)

    if not parse_result.success:
        logger.warning("RDKit parsing failed for SMILES '%s': %s", clean_smiles, parse_result.error_message)
        return MoleculeEvaluationResponse(
            smiles=clean_smiles,
            valid_molecule=False,
            error_detail=parse_result.error_message,
            warnings=pre_check.warnings + parse_result.warnings,
        )

    mol = parse_result.mol

    # -----------------------------------------------------------------------
    # Step 3 — Molecule complexity checks
    # -----------------------------------------------------------------------
    complexity_warnings = check_molecule_complexity(mol, max_atoms=settings.MAX_ATOM_COUNT)
    all_warnings = pre_check.warnings + parse_result.warnings + complexity_warnings

    # -----------------------------------------------------------------------
    # Step 4 — Descriptor computation
    # -----------------------------------------------------------------------
    try:
        raw = compute_descriptors(mol)
    except Exception as exc:
        logger.exception("Descriptor computation failed for SMILES '%s': %s", clean_smiles, exc)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=(
                f"Descriptor computation encountered an unexpected error: {exc}. "
                "Please report this SMILES string to the API maintainer."
            ),
        ) from exc

    # -----------------------------------------------------------------------
    # Step 5 — Lipinski Ro5 evaluation
    # -----------------------------------------------------------------------
    lipinski_result = evaluate_lipinski(
        molecular_weight=raw.molecular_weight,
        logp=raw.logp,
        hbond_donors=raw.hbond_donors,
        hbond_acceptors=raw.hbond_acceptors,
    )
    raw.lipinski_violations = lipinski_result.violations

    # Append Lipinski violation details as informational warnings
    if lipinski_result.violation_details:
        all_warnings.extend(lipinski_result.violation_details)

    # -----------------------------------------------------------------------
    # Step 6 — Medicinal chemistry reasoning
    # -----------------------------------------------------------------------
    analysis_dict = build_analysis(raw, lipinski_result)

    # -----------------------------------------------------------------------
    # Step 7 — Response assembly
    # -----------------------------------------------------------------------
    descriptors = MolecularDescriptors(
        molecular_weight=raw.molecular_weight,
        logp=raw.logp,
        tpsa=raw.tpsa,
        hbond_donors=raw.hbond_donors,
        hbond_acceptors=raw.hbond_acceptors,
        rotatable_bonds=raw.rotatable_bonds,
        lipinski_violations=raw.lipinski_violations,
        heavy_atom_count=raw.heavy_atom_count,
        ring_count=raw.ring_count,
        aromatic_rings=raw.aromatic_rings,
        fraction_csp3=raw.fraction_csp3,
        molar_refractivity=raw.molar_refractivity,
    )

    analysis = MedicinalChemistryAnalysis(**analysis_dict)

    logger.info(
        "Evaluation complete — MW=%.2f LogP=%.2f TPSA=%.1f Ro5=%d violations",
        raw.molecular_weight,
        raw.logp,
        raw.tpsa,
        raw.lipinski_violations,
    )

    return MoleculeEvaluationResponse(
        smiles=clean_smiles,
        canonical_smiles=parse_result.canonical_smiles,
        valid_molecule=True,
        descriptors=descriptors,
        analysis=analysis,
        warnings=all_warnings,
    )


# ---------------------------------------------------------------------------
# GET /examples — curated molecule library for API exploration
# ---------------------------------------------------------------------------

@router.get(
    "/examples",
    summary="Curated SMILES examples for API exploration",
    description=(
        "Returns a library of well-known drug molecules with their SMILES strings, "
        "useful for testing and demonstrating the evaluation engine."
    ),
)
async def get_examples() -> dict:
    """Return a curated set of drug molecules for API testing."""
    return {
        "examples": [
            {
                "name": "Aspirin",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "category": "NSAID / analgesic",
                "notes": "Classic Ro5-compliant, orally bioavailable drug.",
            },
            {
                "name": "Ibuprofen",
                "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
                "category": "NSAID",
                "notes": "Profen scaffold; low MW, good permeability.",
            },
            {
                "name": "Caffeine",
                "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                "category": "Stimulant / xanthine",
                "notes": "Good CNS penetration example.",
            },
            {
                "name": "Atorvastatin (Lipitor)",
                "smiles": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
                "category": "Statin / cholesterol-lowering",
                "notes": "Borderline Ro5 — large but orally bioavailable due to transport.",
            },
            {
                "name": "Paclitaxel (Taxol)",
                "smiles": "O=C(O[C@@H]1C[C@]2(O)C(=O)[C@H](OC(=O)c3ccccc3)[C@@H](O)[C@@H]2[C@@H](OC(=O)[C@@H](NC(=O)c2ccccc2)[C@@H](O)c2ccccc2)[C@H]1OC(=O)C)C",
                "category": "Anticancer / taxane",
                "notes": "Extreme Ro5 violator — highlights limits of Ro5 for natural products.",
            },
            {
                "name": "Penicillin G",
                "smiles": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
                "category": "Antibiotic / beta-lactam",
                "notes": "Good MW but high TPSA — relies on active transport.",
            },
            {
                "name": "Sildenafil (Viagra)",
                "smiles": "CCCC1=NN(C)C(=C1C(=O)NCC)c1cc(S(=O)(=O)N2CCN(C)CC2)ccc1OCC",
                "category": "PDE5 inhibitor",
                "notes": "High MW, moderate Ro5 compliance.",
            },
        ]
    }
