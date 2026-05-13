"""
app/models/schemas.py
=====================
Pydantic v2 data models for all API request and response payloads.

Strict typing and rich field-level documentation ensure that the
auto-generated OpenAPI spec is useful to downstream consumers (front-end
teams, ML pipelines, partner integrations).

Chemistry context is embedded in `description` metadata so that the
interactive /docs UI acts as a lightweight scientific reference.
"""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field, field_validator


# ---------------------------------------------------------------------------
# Request models
# ---------------------------------------------------------------------------


class SMILESRequest(BaseModel):
    """Payload submitted to the descriptor calculation endpoint."""

    smiles: str = Field(
        ...,
        min_length=1,
        max_length=2_000,
        description=(
            "A canonical or non-canonical SMILES string representing the molecule "
            "to be evaluated. Salts and stereochemistry are preserved. "
            "Example: 'CC(=O)OC1=CC=CC=C1C(=O)O' (Aspirin)."
        ),
        examples=["CC(=O)OC1=CC=CC=C1C(=O)O"],
    )

    @field_validator("smiles")
    @classmethod
    def strip_whitespace(cls, v: str) -> str:
        return v.strip()


# ---------------------------------------------------------------------------
# Descriptor sub-model
# ---------------------------------------------------------------------------


class MolecularDescriptors(BaseModel):
    """
    Core physicochemical descriptors derived from the molecular graph.

    Values are computed exclusively from 2-D topology (no 3-D conformer
    generation required), making them suitable for high-throughput virtual
    screening pipelines.
    """

    molecular_weight: float = Field(
        ...,
        ge=0,
        description=(
            "Exact molecular weight in g/mol. "
            "Lipinski's Rule of Five threshold: ≤ 500 Da."
        ),
    )
    logp: float = Field(
        ...,
        description=(
            "Wildman–Crippen partition coefficient (log P). "
            "Quantifies lipophilicity. Lipinski threshold: ≤ 5. "
            "Optimal oral-absorption window: 0–3."
        ),
    )
    tpsa: float = Field(
        ...,
        ge=0,
        description=(
            "Topological Polar Surface Area in Å². "
            "Computed over N, O, S and attached H atoms. "
            "< 90 Å² favours intestinal absorption; < 60 Å² favours BBB penetration."
        ),
    )
    hbond_donors: int = Field(
        ...,
        ge=0,
        description="Number of hydrogen-bond donors (NH + OH groups). Lipinski threshold: ≤ 5.",
    )
    hbond_acceptors: int = Field(
        ...,
        ge=0,
        description=(
            "Number of hydrogen-bond acceptors (N + O atoms). Lipinski threshold: ≤ 10."
        ),
    )
    rotatable_bonds: int = Field(
        ...,
        ge=0,
        description=(
            "Number of rotatable bonds. "
            "Proxy for molecular flexibility and conformational entropy penalty. "
            "Veber's rule: ≤ 10 for good oral bioavailability."
        ),
    )
    lipinski_violations: int = Field(
        ...,
        ge=0,
        le=4,
        description=(
            "Count of Lipinski Ro5 violations (0–4). "
            "≤ 1 violation: drug-like; ≥ 2 violations: poor oral bioavailability predicted."
        ),
    )
    # Extended descriptors (computed but listed for completeness)
    heavy_atom_count: int = Field(..., ge=0, description="Total count of non-hydrogen atoms.")
    ring_count: int = Field(..., ge=0, description="Number of ring systems (SSSR).")
    aromatic_rings: int = Field(..., ge=0, description="Number of aromatic ring systems.")
    fraction_csp3: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description=(
            "Fraction of sp³ carbon atoms. "
            "Higher values correlate with reduced flat-aromatic character, "
            "improved aqueous solubility, and reduced CYP inhibition risk. "
            "Egan / Pfizer: Fsp³ > 0.25 is desirable."
        ),
    )
    molar_refractivity: float = Field(
        ...,
        description=(
            "Molar refractivity (MR). "
            "Encodes molecular volume and polarisability. "
            "Ghose filter: 40 ≤ MR ≤ 130."
        ),
    )


# ---------------------------------------------------------------------------
# Analysis / reasoning sub-model
# ---------------------------------------------------------------------------


class MedicinalChemistryAnalysis(BaseModel):
    """
    Structured medicinal chemistry interpretation of the computed descriptors.

    Each field encodes expert-level reasoning about a pharmacokinetic or
    drug-likeness dimension.  The `optimization_conflicts` list surfaces
    competing physicochemical pressures that must be managed during lead
    optimisation.
    """

    absorption: str = Field(
        ...,
        description=(
            "Predicted oral absorption profile based on MW, TPSA, and rotatable bonds."
        ),
    )
    permeability: str = Field(
        ...,
        description=(
            "Membrane permeability assessment derived from LogP and TPSA. "
            "Relevant to both passive transcellular and paracellular transport."
        ),
    )
    solubility: str = Field(
        ...,
        description="Aqueous solubility prediction based on LogP and TPSA.",
    )
    bbb_penetration: str = Field(
        ...,
        description=(
            "Blood-brain barrier penetration likelihood. "
            "Requires TPSA < 60–90 Å², MW < 450, and moderate LogP."
        ),
    )
    drug_likeness: str = Field(
        ...,
        description="Overall drug-likeness assessment integrating all descriptor dimensions.",
    )
    flexibility_risk: str = Field(
        ...,
        description=(
            "Conformational entropy penalty and binding-pose uncertainty risk "
            "based on rotatable bond count."
        ),
    )
    optimization_conflicts: List[str] = Field(
        default_factory=list,
        description=(
            "Competing physicochemical pressures identified by the reasoning engine. "
            "These represent lead-optimisation tensions that require strategic trade-offs."
        ),
    )
    overall_verdict: str = Field(
        ...,
        description=(
            "Concise synthesis of the molecular profile for medicinal chemistry decision-making."
        ),
    )


# ---------------------------------------------------------------------------
# Top-level response model
# ---------------------------------------------------------------------------


class MoleculeEvaluationResponse(BaseModel):
    """
    Complete evaluation response returned by the /evaluate endpoint.

    When `valid_molecule` is False, only `smiles`, `valid_molecule`, and
    `error_detail` are populated; all other fields are None.
    """

    smiles: str = Field(..., description="The submitted SMILES string (whitespace-stripped).")
    canonical_smiles: Optional[str] = Field(
        None,
        description="RDKit-canonicalised SMILES for unambiguous molecular identity.",
    )
    valid_molecule: bool = Field(
        ...,
        description="True if RDKit successfully parsed the SMILES into a valid molecular graph.",
    )
    error_detail: Optional[str] = Field(
        None,
        description="Human-readable parse or validation error (only present when valid_molecule is False).",
    )
    descriptors: Optional[MolecularDescriptors] = Field(
        None,
        description="Computed physicochemical descriptor profile.",
    )
    analysis: Optional[MedicinalChemistryAnalysis] = Field(
        None,
        description="Medicinal chemistry reasoning derived from the descriptor profile.",
    )
    warnings: List[str] = Field(
        default_factory=list,
        description="Non-fatal chemistry warnings (e.g., unusual valence, disconnected fragments).",
    )

    model_config = {"json_schema_extra": {
        "example": {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "valid_molecule": True,
            "error_detail": None,
            "descriptors": {
                "molecular_weight": 180.16,
                "logp": 1.19,
                "tpsa": 63.6,
                "hbond_donors": 1,
                "hbond_acceptors": 3,
                "rotatable_bonds": 3,
                "lipinski_violations": 0,
                "heavy_atom_count": 13,
                "ring_count": 1,
                "aromatic_rings": 1,
                "fraction_csp3": 0.08,
                "molar_refractivity": 45.26,
            },
            "analysis": {
                "absorption": "Favourable oral absorption predicted …",
                "permeability": "Moderate membrane permeability …",
                "solubility": "Adequate aqueous solubility …",
                "bbb_penetration": "Marginal BBB penetration …",
                "drug_likeness": "Drug-like profile …",
                "flexibility_risk": "Low conformational entropy penalty …",
                "optimization_conflicts": [],
                "overall_verdict": "Aspirin presents a clean drug-like profile …",
            },
            "warnings": [],
        }
    }}
