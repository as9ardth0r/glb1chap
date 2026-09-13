"""glb1chap — pipeline de génération et de criblage de chaperons pharmacologiques
candidats pour la bêta-galactosidase lysosomale humaine (GLB1), déficiente
dans la gangliosidose GM1 et la maladie de Morquio B.

Approche : décoration du scaffold iminosucre DGJ (1-désoxygalactonojirimycine),
un inhibiteur/chaperon de galactosidases déjà co-cristallisé avec GLB1
(PDB 3THD) — même principe pharmacologique que migalastat pour la
Fabry (GLA), appliqué ici à GLB1.
"""

from .properties import MoleculeRecord, compute_batch, compute_descriptors, load_smiles
from .filters import TPPProfile, enrich_and_filter

__all__ = [
    "MoleculeRecord",
    "compute_batch",
    "compute_descriptors",
    "load_smiles",
    "TPPProfile",
    "enrich_and_filter",
]

__version__ = "0.1.0"
