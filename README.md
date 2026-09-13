# glb1chap

> **Génération et criblage in silico de chaperons pharmacologiques candidats pour GLB1 (bêta-galactosidase lysosomale)**

Pipeline open source explorant des dérivés N-substitués du DGJ (1-désoxygalactonojirimycine), un iminosucre déjà connu comme inhibiteur/chaperon de galactosidases et co-cristallisé avec la bêta-galactosidase humaine (PDB [3THD](https://www.rcsb.org/structure/3THD)). L'objectif est d'identifier des candidats respectant un profil cible (TPP) adapté aux chaperons pharmacologiques, dans une optique de recherche sur la gangliosidose GM1 et la maladie de Morquio B — deux maladies de surcharge lysosomale causées par un déficit en GLB1.

L'approche s'inspire directement de celle validée cliniquement pour la Fabry (migalastat, chaperon de l'alpha-galactosidase A/GLA) : une petite molécule polaire portant une amine basique protonable, capable de se fixer dans le site actif d'une enzyme mutante instable et de la stabiliser suffisamment pour restaurer son trafic vers le lysosome.

**Pas de génération libre de l'espace chimique** : le générateur (`generator.py`) se limite à une décoration combinatoire d'un catalogue fini d'une vingtaine de substituants sur l'azote endocyclique du scaffold DGJ — le pharmacophore polyhydroxylé (essentiel à la reconnaissance du site actif) reste toujours intact.

---

## Structure du dépôt

```text
glb1chap/
├── src/glb1chap/
│   ├── properties.py      # Chargement SMILES + calcul MW/LogP/TPSA/HBD/HBA/QED
│   ├── filters.py         # TPP chaperon (amine basique requise, plages MW/LogP/TPSA...), PAINS & Brenk
│   ├── generator.py       # Décoration N-substituée du scaffold DGJ (catalogue fini de R-groups)
│   ├── receptor_prep.py   # Repérage/extraction du ligand co-cristallisé (DGJ) dans un PDB brut
│   ├── docking_prep.py    # Conformères 3D (ETKDGv3 + MMFF94)
│   ├── docking.py         # Docking contre GLB1 (PDBQT, AutoDock Vina)
│   ├── novelty.py         # Vérification PubChem (InChIKey exact)
│   └── hall_of_fame.py    # Persistance des meilleurs candidats entre deux runs
├── scripts/
│   ├── prepare_receptor.py # Télécharge et prépare 3THD pour le docking
│   └── run_pipeline.py     # Orchestration : génération -> filtrage -> nouveauté -> docking -> persistance
├── data/
│   ├── hall_of_fame.json   # Meilleurs candidats conformes au TPP
│   ├── explored.json       # SMILES canoniques déjà testés
│   └── receptor/GLB1/      # Récepteur préparé (PDBQT + boîte Vina), après prepare_receptor.py
├── tests/
└── requirements.txt
```

## Usage

```bash
pip install -r requirements.txt

# 1. Préparer le récepteur (nécessite un accès réseau à files.rcsb.org)
python scripts/prepare_receptor.py

# 2. Lancer le pipeline (génération + filtrage, docking et vérif. PubChem optionnels)
python scripts/run_pipeline.py --dock --check-novelty
```

## Limites connues

- Le catalogue de R-groups est volontairement petit et fixe (pas de génération combinatoire ouverte type mutation atomique) — étendre `generator.R_GROUPS` pour explorer plus large.
- Le SA score est une heuristique simplifiée (pas le SA score officiel RDKit Contrib d'Ertl & Schuffenhauer).
- Le docking (`--dock`) nécessite AutoDock Vina, Meeko et le récepteur déjà préparé ; sans ces prérequis, `docking_score` reste `None` et le tri se fait sur le QED.
- Plusieurs analogues N-substitués du DGJ générés ici (ex. N-butyl-DGJ) sont déjà des composés connus/publiés — `--check-novelty` le signalera via `is_novel=False`, ce qui est une information utile, pas un échec du pipeline.
