"""
glb1chap.docking
==============
Docking moléculaire réel contre le récepteur GLB1 préparé (3THD, voir
scripts/prepare_receptor.py) :
  - préparation du récepteur en PDBQT (fait une fois, en amont),
  - préparation du ligand en PDBQT via Meeko,
  - docking avec AutoDock Vina (bindings Python officiels).
"""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem

TARGET_NAME = "GLB1"


def load_box_config(receptor_dir: str | Path) -> dict | None:
    """Lit la configuration de boîte Vina écrite par prepare_receptor.py
    (fichier <name>_box.txt produit par meeko, format 'key = value')."""
    box_path = Path(receptor_dir) / "glb1_box.txt"
    if not box_path.exists():
        return None
    config: dict[str, float] = {}
    for line in box_path.read_text().splitlines():
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        try:
            config[key.strip()] = float(value.strip())
        except ValueError:
            continue
    if not all(k in config for k in ("center_x", "center_y", "center_z")):
        return None
    return {
        "center": (config["center_x"], config["center_y"], config["center_z"]),
        "box_size": (
            config.get("size_x", 20.0),
            config.get("size_y", 20.0),
            config.get("size_z", 20.0),
        ),
    }


def prepare_ligand_pdbqt(mol: Chem.Mol) -> str | None:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    try:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol)
        pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
        return pdbqt_string if is_ok else None
    except Exception:
        return None


def dock(
    receptor_pdbqt_path: str | Path,
    ligand_pdbqt_string: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
) -> float | None:
    """Lance AutoDock Vina contre le récepteur GLB1 et retourne le meilleur
    score d'affinité prédit (kcal/mol — plus négatif = liaison plus
    favorable). Retourne None en cas d'échec plutôt que de lever une
    exception, pour ne pas interrompre le criblage d'un lot entier."""
    from vina import Vina

    try:
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(str(receptor_pdbqt_path))
        v.set_ligand_from_string(ligand_pdbqt_string)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        scores = v.energies()
        if scores is None or len(scores) == 0:
            return None
        return round(float(scores[0][0]), 2)
    except Exception:
        return None


def dock_smiles(smiles: str, mol_id: str, receptor_dir: str | Path) -> float | None:
    """Enchaîne embedding 3D + préparation ligand + docking pour un SMILES,
    contre le récepteur GLB1 préparé dans `receptor_dir`. Retourne None si
    une étape échoue ou si le récepteur n'est pas préparé (voir
    scripts/prepare_receptor.py)."""
    from .docking_prep import embed_3d

    receptor_dir = Path(receptor_dir)
    receptor_pdbqt = receptor_dir / "glb1.pdbqt"
    if not receptor_pdbqt.exists():
        return None

    box = load_box_config(receptor_dir)
    if box is None:
        return None

    mol = embed_3d(smiles, mol_id=mol_id)
    if mol is None:
        return None

    ligand_pdbqt = prepare_ligand_pdbqt(mol)
    if ligand_pdbqt is None:
        return None

    return dock(receptor_pdbqt, ligand_pdbqt, center=box["center"], box_size=box["box_size"])
