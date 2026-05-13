# MolecularMind — Computational Molecular Evaluation Engine

> A production-grade backend for physicochemical descriptor computation and medicinal chemistry reasoning, built for AI-driven drug discovery workflows.

[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.111-green.svg)](https://fastapi.tiangolo.com/)
[![RDKit](https://img.shields.io/badge/RDKit-2023.9-orange.svg)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

MolecularMind is not a chemistry calculator.

It is a **computational molecular evaluation engine** — a backend system that accepts a SMILES string and returns a richly annotated pharmacokinetic profile that mirrors the kind of analysis performed in real hit-to-lead optimisation campaigns.

The engine computes core physicochemical descriptors using RDKit, evaluates Lipinski's Rule of Five, and applies a multi-dimensional reasoning layer to produce structured interpretations of:

- **Oral absorption** potential
- **Membrane permeability** profile
- **Aqueous solubility** risk
- **Blood-brain barrier** penetration likelihood
- **Conformational flexibility** and entropy penalties
- **Lead-optimisation conflicts** — competing physicochemical pressures

This architecture is designed to serve as the foundation for AI-assisted drug discovery pipelines, ADMET screening tools, and medicinal chemistry decision-support systems.

---

## Architecture

```
backend/
│
├── app/
│   ├── main.py                   # FastAPI app factory, middleware, lifecycle
│   ├── config.py                 # Pydantic BaseSettings — env-driven configuration
│   │
│   ├── routes/
│   │   └── descriptors.py        # POST /evaluate, GET /examples
│   │
│   ├── services/
│   │   ├── rdkit_engine.py       # SMILES parsing & sanitisation (RDKit interface)
│   │   ├── descriptor_calculator.py  # Physicochemical descriptor computation
│   │   ├── medicinal_reasoning.py    # Chemistry reasoning & conflict detection
│   │   └── lipinski.py           # Lipinski Rule of Five evaluator
│   │
│   ├── models/
│   │   └── schemas.py            # Pydantic v2 request/response models
│   │
│   └── utils/
│       └── validators.py         # Pre-parse SMILES validation utilities
│
├── requirements.txt
├── Dockerfile
├── render.yaml
├── .env.example
├── .gitignore
└── README.md
```

### Service Layer Design

```
Request → [validators] → [rdkit_engine] → [descriptor_calculator]
                                                    ↓
                                             [lipinski]
                                                    ↓
                                        [medicinal_reasoning]
                                                    ↓
                                         Pydantic Response Model
```

Each service is a pure-function module with no shared state, making the pipeline fully testable and horizontally scalable.

---

## Descriptors Computed

| Descriptor | Module | Medicinal Chemistry Significance |
|---|---|---|
| Molecular Weight | `Descriptors.ExactMolWt` | Membrane permeability, Ro5 Rule 1 |
| LogP (Crippen) | `Crippen.MolLogP` | Lipophilicity, solubility, permeability |
| TPSA | `rdMolDescriptors.CalcTPSA` | Oral absorption, BBB penetration |
| H-Bond Donors | `CalcNumHBD` | Membrane desolvation penalty |
| H-Bond Acceptors | `CalcNumHBA` | Hydrophilicity, absorption |
| Rotatable Bonds | `CalcNumRotatableBonds` | Conformational entropy, Veber rule |
| Fsp³ | `CalcFractionCSP3` | Saturation, selectivity, solubility |
| Molar Refractivity | `Crippen.MolMR` | Ghose filter, molecular volume |
| Ring Count | `CalcNumRings` | Structural complexity |
| Aromatic Rings | `CalcNumAromaticRings` | PAINS risk, flat-aromatic character |
| Heavy Atom Count | `mol.GetNumHeavyAtoms()` | Molecular size |
| Lipinski Violations | Derived | Drug-likeness assessment |

---

## Scientific Rationale

### Why These Descriptors?

These 12 descriptors represent the **minimum sufficient set** for a first-pass ADME evaluation:

- **Ro5 descriptors** (MW, LogP, HBD, HBA) — Lipinski's Rule of Five remains the most validated filter for oral bioavailability, supported by data from thousands of approved drugs.
- **TPSA** — Ertl's TPSA is superior to calculated log D for absorption and BBB prediction at the virtual screening stage.
- **Rotatable bonds** — Veber et al. (2002) demonstrated that RotBonds and TPSA together outperform MW alone for predicting rat oral bioavailability.
- **Fsp³** — Lovering et al. (2009) showed that increasing molecular saturation (higher Fsp³) correlates with increased clinical success rates and reduced attrition.
- **Molar Refractivity** — The Ghose filter (160 ≤ MW ≤ 480, −0.4 ≤ LogP ≤ 5.6, 40 ≤ MR ≤ 130) complements Ro5 by encoding polarisability.

### Reasoning Engine Philosophy

Real medicinal chemistry is about **managing competing pressures**, not passing/failing filters.

MolecularMind surfaces these tensions explicitly. For example:
- Increasing LogP improves membrane permeability but reduces aqueous solubility and increases CYP inhibition risk.
- Reducing TPSA to improve BBB penetration may remove hydrogen-bond acceptors essential for target engagement.
- Adding sp³ centres (increasing Fsp³) can lower TPSA due to conformational dynamics, creating unexpected solubility effects.

The `optimization_conflicts` field in the response is the most actionable output — it tells a medicinal chemist *where the molecule is hard to optimise*.

---

## Setup & Installation

### Prerequisites

- Python 3.11 or 3.12
- `pip` ≥ 23.0

### Local Development

```bash
# 1. Clone the repository
git clone https://github.com/your-org/molecularmind-backend.git
cd molecularmind-backend

# 2. Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# 4. Configure environment
cp .env.example .env
# Edit .env as needed (RELOAD=true for development)

# 5. Start the development server
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

The API will be live at `http://localhost:8000`.

Interactive documentation: `http://localhost:8000/docs`

### Docker

```bash
# Build
docker build -t molecularmind:latest .

# Run
docker run -p 8000:8000 --env-file .env molecularmind:latest

# Health check
curl http://localhost:8000/health
```

---

## API Reference

### Base URL

```
Production: https://molecularmind-api.onrender.com
Local:      http://localhost:8000
```

### Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/health` | Liveness probe |
| `GET` | `/docs` | Interactive Swagger UI |
| `GET` | `/redoc` | ReDoc documentation |
| `POST` | `/api/v1/evaluate` | Evaluate a molecule from SMILES |
| `GET` | `/api/v1/examples` | Curated SMILES library |

---

### `POST /api/v1/evaluate`

**Request**

```json
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
}
```

**Response — Aspirin**

```json
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "valid_molecule": true,
  "error_detail": null,
  "descriptors": {
    "molecular_weight": 180.0423,
    "logp": 1.3101,
    "tpsa": 63.6,
    "hbond_donors": 1,
    "hbond_acceptors": 3,
    "rotatable_bonds": 3,
    "lipinski_violations": 0,
    "heavy_atom_count": 13,
    "ring_count": 1,
    "aromatic_rings": 1,
    "fraction_csp3": 0.1111,
    "molar_refractivity": 45.2647
  },
  "analysis": {
    "absorption": "The low molecular weight (180.04 Da) provides an excellent foundation for passive intestinal absorption. TPSA of 63.6 Å² is in the optimal oral-absorption window (< 90 Å²).",
    "permeability": "LogP = 1.31 sits in the ideal window for membrane permeability. Adequate lipophilicity for passive transcellular diffusion without excessive non-specific binding or solubility penalties.",
    "solubility": "LogP = 1.31 suggests adequate solubility for most oral formulations.",
    "bbb_penetration": "TPSA = 63.6 Å² allows moderate BBB penetration (60–90 Å² borderline zone). MW = 180.04 Da is within the preferred CNS molecular weight range (< 400 Da). LogP is within the optimal CNS lipophilicity window (1–3).",
    "drug_likeness": "The molecule satisfies all four Lipinski Ro5 criteria with zero violations — an excellent starting point for an orally bioavailable small molecule drug. Fsp³ = 0.11 is in the intermediate range.",
    "flexibility_risk": "Rotatable bonds = 3 — very rigid scaffold. Low conformational entropy penalty upon protein binding.",
    "optimization_conflicts": [],
    "overall_verdict": "Overall drug-likeness: EXCELLENT. MW = 180.0 Da | LogP = 1.31 | TPSA = 63.6 Å² | HBD/HBA = 1/3 | RotBonds = 3 | Fsp³ = 0.11 | Ro5 violations = 0. This compound is a strong hit-to-lead candidate."
  },
  "warnings": []
}
```

**Response — Invalid SMILES**

```json
{
  "smiles": "C(INVALID)",
  "canonical_smiles": null,
  "valid_molecule": false,
  "error_detail": "RDKit was unable to parse the SMILES string. The notation may be malformed.",
  "descriptors": null,
  "analysis": null,
  "warnings": []
}
```

---

### Example API Calls

**cURL**

```bash
# Evaluate aspirin
curl -X POST "http://localhost:8000/api/v1/evaluate" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}'

# Evaluate caffeine
curl -X POST "http://localhost:8000/api/v1/evaluate" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C"}'

# Get example molecules
curl "http://localhost:8000/api/v1/examples"
```

**Python (httpx)**

```python
import httpx

client = httpx.Client(base_url="http://localhost:8000")

response = client.post("/api/v1/evaluate", json={
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
})

data = response.json()
print(f"MW: {data['descriptors']['molecular_weight']}")
print(f"LogP: {data['descriptors']['logp']}")
print(f"Lipinski violations: {data['descriptors']['lipinski_violations']}")
print(f"Verdict: {data['analysis']['overall_verdict']}")
```

**JavaScript (fetch)**

```javascript
const response = await fetch("http://localhost:8000/api/v1/evaluate", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({ smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C" }),
});
const molecule = await response.json();
console.log(molecule.analysis.overall_verdict);
```

---

## Deployment

### Render (Recommended)

1. Push this repository to GitHub.
2. Connect your GitHub repository at [render.com](https://render.com).
3. Select **Blueprint** — Render will auto-detect `render.yaml`.
4. Set `ALLOWED_ORIGINS` in the Render environment variable dashboard.
5. Deploy.

### Railway

```bash
# Install Railway CLI
npm install -g @railway/cli

railway login
railway init
railway up
```

Set environment variables in the Railway dashboard.

### Environment Variables

| Variable | Default | Description |
|---|---|---|
| `HOST` | `0.0.0.0` | Bind address |
| `PORT` | `8000` | Listen port |
| `RELOAD` | `false` | Hot-reload (dev only) |
| `ALLOWED_ORIGINS` | `*` | CORS origins (restrict in prod) |
| `LOG_LEVEL` | `INFO` | Logging verbosity |
| `MAX_ATOM_COUNT` | `500` | Molecule size guard |
| `ENABLE_MORDRED` | `false` | Extended Mordred descriptors |

---

## Future Roadmap

| Feature | Priority | Description |
|---|---|---|
| **Mordred integration** | High | Expand to 1,600+ molecular descriptors via the Mordred library |
| **ESOL solubility** | High | Delaney aqueous solubility model for quantitative solubility prediction |
| **pKa estimation** | High | Ionisation state modelling (Epik / RDKit pKa tools) |
| **Tanimoto similarity** | Medium | Compare submitted molecule against reference compound library |
| **PAINS filter** | Medium | Pan-Assay Interference Compounds structural alert detection |
| **Scaffold decomposition** | Medium | Murcko scaffold extraction and ring system analysis |
| **Synthetic accessibility** | Medium | RDKit SA score for synthesisability estimation |
| **QSAR interface** | Low | scikit-learn model integration for bioactivity prediction |
| **Batch evaluation** | Low | `POST /evaluate/batch` — process arrays of SMILES in one request |
| **3-D conformer generation** | Low | ETKDG conformer + 3-D descriptor computation |
| **Authentication** | Low | API key management for multi-tenant deployment |

---

## Testing

```bash
# Install dev dependencies
pip install pytest pytest-asyncio httpx

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=app --cov-report=html
```

---

## Contributing

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/pka-estimation`.
3. Write tests for all new functionality.
4. Run linting: `ruff check app/ && mypy app/`.
5. Open a pull request with a clear description of the chemistry problem being solved.

---

## Scientific References

1. Lipinski, C. A. et al. (1997). *Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings.* Advanced Drug Delivery Reviews, 23(1–3), 3–25.
2. Veber, D. F. et al. (2002). *Molecular properties that influence the oral bioavailability of drug candidates.* Journal of Medicinal Chemistry, 45(12), 2615–2623.
3. Ertl, P. et al. (2000). *Fast calculation of molecular polar surface area as a sum of fragment-based contributions and its application to the prediction of drug transport properties.* Journal of Medicinal Chemistry, 43(20), 3714–3717.
4. Kelder, J. et al. (1999). *Polar molecular surface as a dominating determinant for oral absorption and brain penetration of drugs.* Pharmaceutical Research, 16(10), 1514–1519.
5. Lovering, F. et al. (2009). *Escape from flatland: Increasing saturation as an approach to improving clinical success.* Journal of Medicinal Chemistry, 52(21), 6752–6756.
6. Wildman, S. A. & Crippen, G. M. (1999). *Prediction of physicochemical parameters by atomic contributions.* Journal of Chemical Information and Computer Sciences, 39(5), 868–873.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

*Built with RDKit, FastAPI, and medicinal chemistry expertise.*
