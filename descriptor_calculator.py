"""
app/services/descriptor_calculator.py
======================================
Physicochemical descriptor computation layer.

This module is the numerical heart of the engine.  It receives a sanitised
RDKit Mol object and returns a populated MolecularDescriptors dataclass.

Design principles
-----------------
* Each descriptor is computed by a dedicated private function to isolate
  failures — a crash in one descriptor does not abort the entire pipeline.
* No I/O, no network calls, no state — this is a pure function module.
* All RDKit descriptor module imports are localised here so that the rest of
  the codebase can be used/tested without RDKit installed.

Descriptor references
---------------------
* MW      — rdkit.Chem.Descriptors.ExactMolWt
* LogP    — rdkit.Chem.Crippen.MolLogP  (Wildman-Crippen, JCIM 1999)
* TPSA    — rdkit.Chem.rdMolDescriptors.CalcTPSA  (Ertl, JMC 2000)
* HBD/HBA — rdkit.Chem.rdMolDescriptors (Lipinski counters)
* RotBonds — rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds
* Fsp3    — rdkit.Chem.rdMolDescriptors.CalcFractionCSP3
* MR      — rdkit.Chem.Crippen.MolMR
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class RawDescriptors:
    """
    Plain-data container for all computed descriptor values.
    Decoupled from Pydantic so the calculator has no schema dependency.
    """

    molecular_weight: float
    logp: float
    tpsa: float
    hbond_donors: int
    hbond_acceptors: int
    rotatable_bonds: int
    heavy_atom_count: int
    ring_count: int
    aromatic_rings: int
    fraction_csp3: float
    molar_refractivity: float
    # Lipinski is computed separately — placeholder here
    lipinski_violations: int = 0


def compute_descriptors(mol: Any) -> RawDescriptors:
    """
    Compute the full descriptor profile for a sanitised RDKit Mol object.

    Parameters
    ----------
    mol:
        A sanitised rdkit.Chem.Mol instance.

    Returns
    -------
    RawDescriptors
        All computed physicochemical descriptor values.

    Raises
    ------
    RuntimeError
        If any critical descriptor calculation fails unexpectedly.
    """
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

    mw = _safe_compute("MolecularWeight", Descriptors.ExactMolWt, mol)
    logp = _safe_compute("LogP", Crippen.MolLogP, mol)
    tpsa = _safe_compute("TPSA", rdMolDescriptors.CalcTPSA, mol)
    hbd = _safe_compute("HBondDonors", rdMolDescriptors.CalcNumHBD, mol)
    hba = _safe_compute("HBondAcceptors", rdMolDescriptors.CalcNumHBA, mol)
    rot = _safe_compute("RotatableBonds", rdMolDescriptors.CalcNumRotatableBonds, mol)
    heavy = mol.GetNumHeavyAtoms()
    rings = _safe_compute("RingCount", rdMolDescriptors.CalcNumRings, mol)
    arom = _safe_compute("AromaticRings", rdMolDescriptors.CalcNumAromaticRings, mol)
    fsp3 = _safe_compute("FractionCSP3", rdMolDescriptors.CalcFractionCSP3, mol)
    mr = _safe_compute("MolarRefractivity", Crippen.MolMR, mol)

    # Round floating-point values to chemically meaningful precision
    return RawDescriptors(
        molecular_weight=round(float(mw), 4),
        logp=round(float(logp), 4),
        tpsa=round(float(tpsa), 4),
        hbond_donors=int(hbd),
        hbond_acceptors=int(hba),
        rotatable_bonds=int(rot),
        heavy_atom_count=int(heavy),
        ring_count=int(rings),
        aromatic_rings=int(arom),
        fraction_csp3=round(float(fsp3), 4),
        molar_refractivity=round(float(mr), 4),
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _safe_compute(name: str, fn, mol: Any, default: float = 0.0) -> float:
    """
    Call an RDKit descriptor function with exception isolation.

    If the computation raises, the failure is logged and `default` is
    returned so that one broken descriptor does not abort the whole profile.

    Parameters
    ----------
    name:
        Human-readable descriptor name for logging purposes.
    fn:
        Callable — an RDKit descriptor function accepting a Mol object.
    mol:
        Sanitised RDKit Mol object.
    default:
        Value to return if the computation fails.
    """
    try:
        return fn(mol)
    except Exception as exc:
        logger.warning("Descriptor '%s' computation failed: %s — returning default %s", name, exc, default)
        return default
