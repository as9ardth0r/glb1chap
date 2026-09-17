# glb1chap

> **In silico generation and screening of candidate pharmacological chaperones for GLB1 (lysosomal beta-galactosidase)**

Open-source pipeline exploring N-substituted derivatives of DGJ
(1-deoxygalactonojirimycin), an iminosugar already known as a
galactosidase inhibitor/chaperone and co-crystallized with human
beta-galactosidase (PDB [3THD](https://www.rcsb.org/structure/3THD)).
The goal is to identify candidates meeting a target product profile
(TPP) suited to pharmacological chaperones, in the context of research
on GM1 gangliosidosis and Morquio B disease — two lysosomal storage
disorders caused by GLB1 deficiency.

The approach directly mirrors the one clinically validated for Fabry
disease (migalastat, a chaperone for alpha-galactosidase A/GLA): a
small polar molecule carrying a protonatable basic amine, able to bind
the active site of an unstable mutant enzyme and stabilize it enough
to restore its trafficking to the lysosome.

**Two-stage generation, indefinitely renewable:**
1. **Catalog** (`generator.R_GROUPS`, ~40 substituents) — finite,
   exhausted after the first complete run.
2. **Atomic mutation** (`generator.mutate_fragment`) — once the catalog
   is exhausted, the pipeline mutates the R fragment of the best
   candidates from the hall of fame (adding/removing/permuting a
   halogen, adding/removing a methyl group). This mutation always
   operates **on the isolated fragment, never on the assembled
   molecule**: the DGJ's polyhydroxylated core (the pharmacophore
   recognized by the active site) can structurally never be touched,
   even after hundreds of generations. This mechanism is what lets the
   pipeline run indefinitely (see `.github/workflows/daily-run.yml`,
   an automatic daily run).

---

## Repository structure

```text
glb1chap/
├── src/glb1chap/
│   ├── properties.py      # SMILES loading + MW/LogP/TPSA/HBD/HBA/QED calculation
│   ├── filters.py         # Chaperone TPP (required basic amine, MW/LogP/TPSA ranges...), PAINS & Brenk
│   ├── generator.py       # R-group catalog + atomic fragment mutation (open-ended space)
│   ├── evolve.py          # Combines catalog and mutation to produce a run's candidates
│   ├── receptor_prep.py   # Locating/extracting the co-crystallized ligand (DGJ) from a raw PDB
│   ├── docking_prep.py    # 3D conformers (ETKDGv3 + MMFF94) + SDF export for the 3D viewer
│   ├── docking.py         # Docking against GLB1 (PDBQT, AutoDock Vina)
│   ├── novelty.py         # PubChem check (exact InChIKey)
│   ├── hall_of_fame.py    # Persistence + parent selection for mutation
│   └── export.py          # Generates site/data/molecules.json and conformers.json
├── scripts/
│   ├── prepare_receptor.py # Downloads and prepares 3THD for docking
│   └── run_pipeline.py     # Full orchestration (see below)
├── site/                   # Static dashboard (2D structures + 3D viewer per molecule)
├── data/
│   ├── hall_of_fame.json   # Best TPP-compliant candidates (also the parent pool for mutation)
│   ├── explored.json       # Canonical SMILES already tested, across all mechanisms
│   └── receptor/GLB1/      # Prepared receptor (PDBQT + Vina box), after prepare_receptor.py
├── .github/workflows/
│   ├── daily-run.yml       # Automatic daily run (cron): generates + commits results
│   └── deploy.yml          # Publishes site/ to GitHub Pages on every committed change
├── tests/
└── requirements.txt
```

## Usage

```bash
pip install -r requirements.txt

# 1. Prepare the receptor (requires network access to files.rcsb.org) — one time only
python scripts/prepare_receptor.py

# 2. Run the pipeline
python scripts/run_pipeline.py --dock --check-novelty --n-mutants 30
```

`--n-mutants` controls how many candidates are produced by mutation on
each run (0 to only test the catalog's remainder). Once the catalog is
exhausted, this mutation is what feeds the pipeline indefinitely.

## Automation

- **`daily-run.yml`** runs every night (cron `0 3 * * *`), executes
  `run_pipeline.py --check-novelty` (docking is skipped by default in
  the automated run — see the comment in the workflow), and commits
  `data/*.json` + `site/data/*.json` if there's anything new.
- **`deploy.yml`** triggers on every push touching `site/` and
  republishes the dashboard to GitHub Pages — it doesn't regenerate
  anything itself, it publishes the state already committed by
  `daily-run.yml` (or a manual run).
- To enable: Settings → Pages → Source → *GitHub Actions*, on the
  GitHub repository.

## Dashboard

Sortable cards (fitness / docking score / QED / MW), a "confirmed
novelty" filter, and a **View in 3D** button per molecule that loads
`site/data/conformers.json` on demand and displays the conformer via
[3Dmol.js](https://3dmol.csb.pitt.edu/). The 🧬 icon on a card
indicates a molecule produced by mutation (with its parent shown in a
tooltip) rather than from the starting catalog.

```bash
cd site && python -m http.server 8000   # then open http://localhost:8000
```
(opening `index.html` directly via `file://` blocks the JSON `fetch`
calls in most browsers)

## Known limitations

- The SA score is a simplified heuristic (not the official RDKit
  Contrib SA score by Ertl & Schuffenhauer).
- Docking (`--dock`) requires AutoDock Vina, Meeko, and an already-
  prepared receptor; without these prerequisites, `docking_score`
  stays `None` and sorting falls back to QED. Not enabled in the
  automated daily run by default.
- Atomic mutation is deliberately simple (halogens + methyls) — no
  ring formation, no scaffold change, no rearrangement.
- Several catalog analogs (e.g. N-butyl-DGJ) are already known/
  published compounds — `--check-novelty` flags this via
  `is_novel=False`, which is useful information, not a pipeline
  failure.
