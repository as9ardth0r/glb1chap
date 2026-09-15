# glb1chap

> **Génération et criblage in silico de chaperons pharmacologiques candidats pour GLB1 (bêta-galactosidase lysosomale)**

Pipeline open source explorant des dérivés N-substitués du DGJ (1-désoxygalactonojirimycine), un iminosucre déjà connu comme inhibiteur/chaperon de galactosidases et co-cristallisé avec la bêta-galactosidase humaine (PDB [3THD](https://www.rcsb.org/structure/3THD)). L'objectif est d'identifier des candidats respectant un profil cible (TPP) adapté aux chaperons pharmacologiques, dans une optique de recherche sur la gangliosidose GM1 et la maladie de Morquio B — deux maladies de surcharge lysosomale causées par un déficit en GLB1.

L'approche s'inspire directement de celle validée cliniquement pour la Fabry (migalastat, chaperon de l'alpha-galactosidase A/GLA) : une petite molécule polaire portant une amine basique protonable, capable de se fixer dans le site actif d'une enzyme mutante instable et de la stabiliser suffisamment pour restaurer son trafic vers le lysosome.

**Génération en deux temps, indéfiniment renouvelable :**
1. **Catalogue** (`generator.R_GROUPS`, ~40 substituants) — fini, épuisé après le premier run complet.
2. **Mutation atomique** (`generator.mutate_fragment`) — une fois le catalogue épuisé, le pipeline mute le fragment R des meilleurs candidats du hall of fame (ajout/retrait/permutation d'un halogène, ajout/retrait d'un méthyle). Cette mutation opère **toujours sur le fragment isolé, jamais sur la molécule assemblée** : le noyau polyhydroxylé du DGJ (le pharmacophore reconnu par le site actif) ne peut structurellement jamais être touché, même après des centaines de générations. C'est ce mécanisme qui permet au pipeline de tourner indéfiniment (voir `.github/workflows/daily-run.yml`, exécution quotidienne automatique).

---

## Structure du dépôt

```text
glb1chap/
├── src/glb1chap/
│   ├── properties.py      # Chargement SMILES + calcul MW/LogP/TPSA/HBD/HBA/QED
│   ├── filters.py         # TPP chaperon (amine basique requise, plages MW/LogP/TPSA...), PAINS & Brenk
│   ├── generator.py       # Catalogue de R-groups + mutation atomique du fragment (espace ouvert)
│   ├── evolve.py          # Combine catalogue et mutation pour produire les candidats d'un run
│   ├── receptor_prep.py   # Repérage/extraction du ligand co-cristallisé (DGJ) dans un PDB brut
│   ├── docking_prep.py    # Conformères 3D (ETKDGv3 + MMFF94) + export SDF pour le viewer 3D
│   ├── docking.py         # Docking contre GLB1 (PDBQT, AutoDock Vina)
│   ├── novelty.py         # Vérification PubChem (InChIKey exact)
│   ├── hall_of_fame.py    # Persistance + sélection des parents pour la mutation
│   └── export.py          # Génère site/data/molecules.json et conformers.json
├── scripts/
│   ├── prepare_receptor.py # Télécharge et prépare 3THD pour le docking
│   └── run_pipeline.py     # Orchestration complète (voir ci-dessous)
├── site/                   # Dashboard statique (structures 2D + viewer 3D par molécule)
├── data/
│   ├── hall_of_fame.json   # Meilleurs candidats conformes au TPP (sert aussi de pool de parents)
│   ├── explored.json       # SMILES canoniques déjà testés, tous mécanismes confondus
│   └── receptor/GLB1/      # Récepteur préparé (PDBQT + boîte Vina), après prepare_receptor.py
├── .github/workflows/
│   ├── daily-run.yml       # Run quotidien automatique (cron) : génère + commite les résultats
│   └── deploy.yml          # Publie site/ sur GitHub Pages à chaque changement committé
├── tests/
└── requirements.txt
```

## Usage

```bash
pip install -r requirements.txt

# 1. Préparer le récepteur (nécessite un accès réseau à files.rcsb.org) — une seule fois
python scripts/prepare_receptor.py

# 2. Lancer le pipeline
python scripts/run_pipeline.py --dock --check-novelty --n-mutants 30
```

`--n-mutants` contrôle combien de candidats sont produits par mutation à chaque run (0 pour ne tester que le reliquat du catalogue). Une fois le catalogue épuisé, c'est cette mutation qui alimente le pipeline indéfiniment.

## Automatisation

- **`daily-run.yml`** tourne chaque nuit (cron `0 3 * * *`), exécute `run_pipeline.py --check-novelty` (le docking est omis par défaut du run automatique — voir commentaire dans le workflow), et commite `data/*.json` + `site/data/*.json` s'il y a du nouveau.
- **`deploy.yml`** se déclenche à chaque push touchant `site/` et republie le dashboard sur GitHub Pages — il ne régénère rien lui-même, il publie l'état déjà committé par `daily-run.yml` (ou par un run manuel).
- Activer : Settings → Pages → Source → *GitHub Actions*, sur le dépôt GitHub.

## Dashboard

Cartes triables (fitness / score de docking / QED / MW), filtre "nouveauté confirmée", et un bouton **Voir en 3D** par molécule qui charge à la demande `site/data/conformers.json` et affiche le conformère via [3Dmol.js](https://3dmol.csb.pitt.edu/). L'icône 🧬 sur une carte indique une molécule issue de la mutation (avec son parent en info-bulle) plutôt que du catalogue de départ.

```bash
cd site && python -m http.server 8000   # puis ouvrir http://localhost:8000
```
(l'ouverture directe de `index.html` en `file://` bloque le `fetch` des JSON dans la plupart des navigateurs)

## Limites connues

- Le SA score est une heuristique simplifiée (pas le SA score officiel RDKit Contrib d'Ertl & Schuffenhauer).
- Le docking (`--dock`) nécessite AutoDock Vina, Meeko et le récepteur déjà préparé ; sans ces prérequis, `docking_score` reste `None` et le tri se fait sur le QED. Non activé dans le run quotidien automatique par défaut.
- La mutation atomique est volontairement simple (halogènes + méthyles) — elle ne fait ni cyclisation, ni changement de squelette, ni réarrangement.
- Plusieurs analogues du catalogue (ex. N-butyl-DGJ) sont déjà des composés connus/publiés — `--check-novelty` le signale via `is_novel=False`, ce qui est une information utile, pas un échec du pipeline.
