"""
app/services/rdkit_engine.py
============================
Low-level RDKit interface — molecular parsing, sanitisation, and object
lifecycle management.

This module is the only place in the codebase that directly imports from
rdkit.Chem.  All higher-level services receive a validated Mol object from
here, which enforces a clean dependency boundary and makes the engine
testable without RDKit (via mocking).

Chemistry note
--------------
RDKit's Chem.MolFromSmiles() performs:
  1. SMILES tokenisation and graph construction
  2. Aromaticity perception (default: Hückel model)
  3. Sanitisation (valence checking, ring-info computation)

We catch sanitisation failures explicitly because they carry diagnostic
information useful for the error response (e.g., "Explicit valence for atom
# 3, N, is greater than permitted").
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class ParseResult:
    """Outcome of an RDKit parse attempt."""

    success: bool
    mol: Optional[Any] = None          # rdkit.Chem.Mol — typed as Any to avoid import-time dep
    canonical_smiles: Optional[str] = None
    error_message: Optional[str] = None
    warnings: List[str] = field(default_factory=list)


def parse_smiles(smiles: str) -> ParseResult:
    """
    Parse a SMILES string into a sanitised RDKit Mol object.

    Sanitisation is run in two stages so that informative error messages
    can be extracted even when the molecule is invalid:

    Stage 1 — Parse without sanitisation to get the raw graph.
    Stage 2 — Run sanitisation steps explicitly, catching any failure.

    Parameters
    ----------
    smiles:
        A pre-validated SMILES string (whitespace already stripped).

    Returns
    -------
    ParseResult
        Contains the Mol object on success or a descriptive error on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import SanitizeMol, SanitizeFlags
    except ImportError as exc:
        logger.critical("RDKit is not available: %s", exc)
        return ParseResult(
            success=False,
            error_message="RDKit is not installed on this server. Contact the administrator.",
        )

    # -----------------------------------------------------------------------
    # Stage 1: parse with sanitisation disabled to get the raw graph
    # -----------------------------------------------------------------------
    mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_raw is None:
        logger.warning("RDKit could not tokenise SMILES: %s", smiles)
        return ParseResult(
            success=False,
            error_message=(
                "RDKit was unable to parse the SMILES string. "
                "The notation may be malformed. Please verify atom symbols, "
                "bond types, ring-closure digits, and bracket atoms."
            ),
        )

    # -----------------------------------------------------------------------
    # Stage 2: explicit sanitisation — catches valence errors, etc.
    # -----------------------------------------------------------------------
    try:
        SanitizeMol(mol_raw, catchErrors=False)
    except Exception as sanitise_err:  # rdkit raises Chem.AtomValenceException, etc.
        err_str = str(sanitise_err)
        logger.warning("Sanitisation failed for SMILES '%s': %s", smiles, err_str)

        # Attempt to give a more user-friendly message for the most common errors
        friendly = _translate_sanitisation_error(err_str)
        return ParseResult(
            success=False,
            error_message=friendly,
        )

    # -----------------------------------------------------------------------
    # Canonicalise for reproducibility
    # -----------------------------------------------------------------------
    canonical = Chem.MolToSmiles(mol_raw, canonical=True)

    logger.debug("Successfully parsed SMILES → canonical: %s", canonical)
    return ParseResult(
        success=True,
        mol=mol_raw,
        canonical_smiles=canonical,
    )


def _translate_sanitisation_error(rdkit_error: str) -> str:
    """
    Convert a raw RDKit error message into a medicinal-chemistry-friendly explanation.

    RDKit error strings can be cryptic for users unfamiliar with cheminformatics.
    This mapping improves developer experience and error transparency.
    """
    error_lower = rdkit_error.lower()

    if "valence" in error_lower:
        return (
            f"Valence error: {rdkit_error}. "
            "One or more atoms have an incorrect number of bonds for their element type. "
            "Common causes: over-specified nitrogen (N with 5 bonds without explicit charge), "
            "or missing hydrogen counts on bracket atoms."
        )
    if "ring" in error_lower and "closure" in error_lower:
        return (
            f"Ring-closure error: {rdkit_error}. "
            "A ring-closure digit was opened but never closed, or vice versa."
        )
    if "aromaticity" in error_lower:
        return (
            f"Aromaticity perception failed: {rdkit_error}. "
            "The aromatic system described by lowercase atoms does not satisfy Hückel's rule. "
            "Use uppercase atoms with explicit double bonds for non-aromatic conjugated systems."
        )
    if "atom" in error_lower and "not" in error_lower:
        return (
            f"Unrecognised atom type: {rdkit_error}. "
            "Check that all element symbols are valid (RDKit is case-sensitive for "
            "two-letter elements in brackets, e.g., [Fe], not [fe])."
        )

    # Fallback — return the raw RDKit message wrapped with context
    return (
        f"Molecular graph sanitisation failed: {rdkit_error}. "
        "Please review the SMILES notation for structural consistency."
    )
