#!/usr/bin/env python3
"""Orchestration du pipeline GLB1 : génère les analogues N-substitués du
DGJ (catalogue fini, voir generator.R_GROUPS), calcule leurs descripteurs,
les filtre contre le TPP chaperon, vérifie leur nouveauté (PubChem), les
docke contre GLB1 (3THD) si le récepteur est préparé, et met à jour le
hall of fame persistant.

Contrairement à un espace combinatoire ouvert, le catalogue de R-groups
est petit et fini : ce script peut simplement le parcourir en entier à
chaque run plutôt que d'échantillonner — pas besoin d'une boucle évolutive.

Usage :
    python scripts/run_pipeline.py [--dock] [--check-novelty]
"""
from __future__ import annotations

import argparse
import sys
from datetime import date, timezone, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from glb1chap import export, hall_of_fame, novelty  # noqa: E402
from glb1chap.filters import enrich_and_filter  # noqa: E402
from glb1chap.generator import generate_all, mol_id_for  # noqa: E402
from glb1chap.properties import compute_descriptors  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RECEPTOR_DIR = DATA_DIR / "receptor" / "GLB1"
SITE_DATA_PATH = REPO_ROOT / "site" / "data" / "molecules.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dock", action="store_true",
                         help="docke chaque candidat conforme au TPP contre GLB1 (3THD) — "
                              "nécessite d'avoir lancé scripts/prepare_receptor.py au préalable")
    parser.add_argument("--check-novelty", action="store_true",
                         help="interroge PubChem par InChIKey pour chaque candidat conforme")
    args = parser.parse_args()

    today = date.today().isoformat()
    hall_path = DATA_DIR / "hall_of_fame.json"
    explored_path = DATA_DIR / "explored.json"

    explored = hall_of_fame.load_explored(explored_path)
    hall = hall_of_fame.load_hall_of_fame(hall_path)

    passing = []
    for smiles, recipe in generate_all():
        mol = compute_descriptors(smiles, mol_id_for(smiles, recipe, today))
        if mol.canonical_smiles in explored:
            continue
        explored.add(mol.canonical_smiles)
        mol.recipe = {"scaffold": recipe.scaffold, "r_group": recipe.r_group}
        mol = enrich_and_filter(mol)
        print(f"[run_pipeline] {recipe.r_group:20s} tpp_pass={mol.tpp_pass}  {mol.notes[:70]}")
        if not mol.tpp_pass:
            continue

        if args.check_novelty:
            inchikey = novelty.inchikey_for(mol.canonical_smiles)
            if inchikey:
                hit = novelty.check_pubchem(inchikey)
                mol.is_novel = hit is None
                if hit:
                    mol.pubchem_cid = hit.get("cid")
                    mol.chemical_name = hit.get("iupac_name")

        if args.dock:
            from glb1chap.docking import dock_smiles
            mol.docking_score = dock_smiles(mol.canonical_smiles, mol.id, RECEPTOR_DIR)
            if mol.docking_score is None:
                print(f"[run_pipeline] docking indisponible pour {recipe.r_group} "
                      f"(récepteur non préparé ? voir scripts/prepare_receptor.py)")

        passing.append(mol)

    hall = hall_of_fame.merge_into_hall_of_fame(hall, passing, today)
    hall_of_fame.save_hall_of_fame(hall, hall_path)
    hall_of_fame.save_explored(explored, explored_path)

    payload = export.build_site_payload(hall)
    export.write_json(payload, SITE_DATA_PATH)

    print(f"\n[run_pipeline] {len(passing)} nouveau(x) candidat(s) conforme(s) au TPP "
          f"sur ce run ; hall of fame : {len(hall)} molécule(s) au total.")
    print(f"[run_pipeline] dashboard mis à jour -> {SITE_DATA_PATH}")


if __name__ == "__main__":
    main()
