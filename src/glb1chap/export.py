"""
glb1chap.export
=============
Assemble les données du hall of fame dans le schéma JSON consommé par le
dashboard statique (site/), avec une dépiction 2D SVG par molécule
(générée côté serveur avec RDKit — pas de dépendance JS de dessin chimique
dans le navigateur).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .properties import MoleculeRecord


def mol_to_svg(smiles: str, width: int = 260, height: int = 200) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def build_site_payload(records: list[MoleculeRecord], target: str = "GLB1") -> dict:
    molecules = []
    for r in records:
        entry = r.to_dict()
        entry["svg"] = mol_to_svg(r.canonical_smiles)
        molecules.append(entry)
    return {
        "target": target,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "count": len(molecules),
        "molecules": molecules,
    }


def write_json(payload: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
