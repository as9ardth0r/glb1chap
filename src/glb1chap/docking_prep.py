"""
glb1chap.docking_prep
====================
Préparation des ligands pour le docking moléculaire : génération de
conformères 3D (ETKDGv3), optimisation MMFF94, export SDF/PDB.

Les candidats iminosucres portent plusieurs stéréocentres fixes (hérités
du scaffold DGJ) : l'embedding doit préserver la stéréochimie du SMILES
d'entrée, pas la deviner.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def embed_3d(smiles: str, mol_id: str = "mol", seed: int = 42) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
        if cid < 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    mol.SetProp("_Name", mol_id)
    return mol


def export_sdf(mol: Chem.Mol, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def export_pdb(mol: Chem.Mol, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(path))
