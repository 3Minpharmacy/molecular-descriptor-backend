"""
app/utils/validators.py
=======================
Input validation utilities for the molecular evaluation engine.

Centralising validation logic here keeps the route handlers and service
layer clean and ensures that every code path that processes a SMILES string
goes through the same safety checks.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import List, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Allowed character set for SMILES strings
# SMILES may contain: atoms, bonds, branches, ring-closure digits,
# stereochemistry descriptors, charges, isotope labels, and whitespace
# (which is stripped before validation).
# ---------------------------------------------------------------------------
_SMILES_ALLOWED = re.compile(
    r"^[A-Za-z0-9@+\-\[\]\(\)\.\#\$\=\/\\%:*,\s]+$"
)

# Characters that hint at common user mistakes (e.g., InChI prefix)
_INCHI_PATTERN = re.compile(r"^InChI=", re.IGNORECASE)
_SMARTS_PATTERN = re.compile(r"[$&,;~!]")  # SMARTS-specific operators


@dataclass
class ValidationResult:
    """Outcome of a pre-parse SMILES validation check."""

    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    sanitised_smiles: Optional[str] = None


def validate_smiles_string(smiles: str, max_atom_hint: int = 500) -> ValidationResult:
    """
    Perform lightweight, RDKit-independent pre-validation of a SMILES string.

    This function is intentionally fast and does NOT attempt to parse the
    molecular graph.  It catches the most common user-input errors before
    handing off to the heavier RDKit parsing stage.

    Parameters
    ----------
    smiles:
        Raw SMILES string submitted by the user.
    max_atom_hint:
        Soft upper bound on expected atom count derived from SMILES length.
        A very long SMILES string is a proxy for a pathologically large molecule.

    Returns
    -------
    ValidationResult
        Structured outcome with error and warning messages.
    """
    result = ValidationResult(is_valid=True)
    cleaned = smiles.strip()

    # --- Emptiness -----------------------------------------------------------
    if not cleaned:
        result.is_valid = False
        result.errors.append("SMILES string is empty.")
        return result

    result.sanitised_smiles = cleaned

    # --- InChI submitted by mistake ------------------------------------------
    if _INCHI_PATTERN.match(cleaned):
        result.is_valid = False
        result.errors.append(
            "An InChI string was submitted. This API accepts SMILES notation only. "
            "Please convert your InChI to SMILES before submitting."
        )
        return result

    # --- SMARTS operators (not valid molecule SMILES) -------------------------
    if _SMARTS_PATTERN.search(cleaned):
        result.warnings.append(
            "The SMILES string contains characters associated with SMARTS query notation "
            "($, &, ;, ~, !). If this is a query pattern rather than a molecule, "
            "results may be unreliable."
        )

    # --- Character set -------------------------------------------------------
    if not _SMILES_ALLOWED.match(cleaned):
        result.is_valid = False
        illegal = sorted(
            {c for c in cleaned if not re.match(r"[A-Za-z0-9@+\-\[\]\(\)\.\#\$\=\/\\%:*,\s]", c)}
        )
        result.errors.append(
            f"SMILES string contains illegal characters: {illegal}. "
            "Only standard SMILES notation is accepted."
        )
        return result

    # --- Length heuristic (molecule size guard) --------------------------------
    if len(cleaned) > max_atom_hint * 4:
        result.warnings.append(
            f"SMILES string is unusually long ({len(cleaned)} characters). "
            "Descriptor calculation may be slow for very large molecules."
        )

    # --- Bracket balance ------------------------------------------------------
    if cleaned.count("[") != cleaned.count("]"):
        result.is_valid = False
        result.errors.append(
            "Unbalanced square brackets detected. "
            "Every opening '[' must have a matching ']'."
        )
        return result

    if cleaned.count("(") != cleaned.count(")"):
        result.is_valid = False
        result.errors.append(
            "Unbalanced parentheses detected. "
            "Every opening '(' must have a matching ')'."
        )
        return result

    # --- Trivial single-atom guard (not necessarily invalid, just warn) -------
    if len(cleaned) == 1 and cleaned.isalpha():
        result.warnings.append(
            "Single-atom SMILES submitted. Descriptor values will be degenerate."
        )

    logger.debug("Pre-validation passed for SMILES: %s", cleaned)
    return result


def check_molecule_complexity(mol, max_atoms: int = 500) -> List[str]:
    """
    Post-parse complexity checks on an RDKit Mol object.

    Returns a list of warning strings (empty if the molecule is unremarkable).

    Parameters
    ----------
    mol:
        A validated RDKit Mol object.
    max_atoms:
        Hard limit on atom count beyond which results are flagged.
    """
    from rdkit.Chem import rdMolDescriptors  # lazy import — rdkit must be available

    warnings: List[str] = []
    atom_count = mol.GetNumAtoms()

    if atom_count > max_atoms:
        warnings.append(
            f"Molecule contains {atom_count} heavy atoms, exceeding the recommended "
            f"limit of {max_atoms}. Descriptor accuracy is not guaranteed for "
            "macromolecules or polymers."
        )

    # Disconnected fragments (mixture / salt)
    from rdkit.Chem import rdmolops
    frags = rdmolops.GetMolFrags(mol)
    if len(frags) > 1:
        warnings.append(
            f"Molecule contains {len(frags)} disconnected fragments (mixture or salt). "
            "Descriptors are computed on the complete graph. Consider desalting "
            "before evaluation for a single-component analysis."
        )

    return warnings
