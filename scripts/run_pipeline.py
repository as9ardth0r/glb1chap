#!/usr/bin/env python3
"""Orchestration du pipeline GLB1, pensée pour un run quotidien récurrent
(voir .github/workflows/daily-run.yml) :

  1. Génère les candidats non encore testés du catalogue fini
     (generator.R_GROUPS) — s'épuise après le premier run complet.
  2. Génère des candidats supplémentaires par mutation atomique des
     meilleurs parents du hall of fame (--n-mutants) — espace ouvert,
     c'est ce qui permet au pipeline de continuer à produire du nouveau
     indéfiniment une fois le catalogue épuisé.
  3. Calcule les descripteurs, filtre contre le TPP chaperon.
  4. Vérifie la nouveauté (PubChem) et docke contre GLB1 (3THD) si demandé.
  5. Met à jour le hall of fame persistant et exporte le dashboard
     (site/data/molecules.json + site/data/conformers.json pour la vue 3D).

Usage :
    python scripts/run_pipeline.py [--dock] [--check-novelty] [--n-mutants 20]
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from glb1chap import evolve, export, hall_of_fame, novelty  # noqa: E402
from glb1chap.filters import enrich_and_filter  # noqa: E402
from glb1chap.generator import mol_id_for  # noqa: E402
from glb1chap.properties import compute_descriptors  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RECEPTOR_DIR = DATA_DIR / "receptor" / "GLB1"
SITE_DATA_DIR = REPO_ROOT / "site" / "data"


def process_candidates(candidates, explored, args, today):
    """Calcule descripteurs + filtre TPP + (optionnel) nouveauté/docking
    pour une liste de (smiles, recipe), en ignorant ceux déjà explorés.
    Retourne les MoleculeRecord conformes au TPP."""
    passing = []
    for smiles, recipe in candidates:
        mol = compute_descriptors(smiles, mol_id_for(smiles, recipe, today))
        if mol.canonical_smiles in explored:
            continue
        explored.add(mol.canonical_smiles)
        mol.recipe = {
            "scaffold": recipe.scaffold,
            "r_group": recipe.r_group,
            "r_group_smiles": recipe.r_group_smiles,
            "parent_id": recipe.parent_id,
        }
        mol = enrich_and_filter(mol)
        print(f"[run_pipeline] {recipe.r_group:24s} tpp_pass={mol.tpp_pass}  {mol.notes[:60]}")
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
                      f"(récepteur préparé ? voir scripts/prepare_receptor.py)")

        passing.append(mol)
    return passing


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dock", action="store_true",
                         help="docke chaque candidat conforme au TPP contre GLB1 (3THD)")
    parser.add_argument("--check-novelty", action="store_true",
                         help="interroge PubChem par InChIKey pour chaque candidat conforme")
    parser.add_argument("--n-mutants", type=int, default=20,
                         help="nombre de candidats à produire par mutation atomique (0 pour désactiver)")
    parser.add_argument("--n-parents", type=int, default=5,
                         help="nombre de meilleurs parents utilisés comme base de mutation")
    args = parser.parse_args()

    today = date.today().isoformat()
    hall_path = DATA_DIR / "hall_of_fame.json"
    explored_path = DATA_DIR / "explored.json"

    explored = hall_of_fame.load_explored(explored_path)
    hall = hall_of_fame.load_hall_of_fame(hall_path)

    catalog_candidates = evolve.candidates_from_catalog(explored)
    print(f"[run_pipeline] {len(catalog_candidates)} candidat(s) du catalogue à évaluer")
    passing = process_candidates(catalog_candidates, explored, args, today)

    if args.n_mutants > 0:
        mutant_candidates = evolve.candidates_from_mutation(
            hall, explored, n_mutants=args.n_mutants, n_parents=args.n_parents,
        )
        print(f"[run_pipeline] {len(mutant_candidates)} candidat(s) muté(s) à évaluer "
              f"(parents pris parmi les {args.n_parents} meilleurs du hall of fame)")
        passing += process_candidates(mutant_candidates, explored, args, today)
    else:
        mutant_candidates = []

    hall = hall_of_fame.merge_into_hall_of_fame(hall, passing, today)
    hall_of_fame.save_hall_of_fame(hall, hall_path)
    hall_of_fame.save_explored(explored, explored_path)

    payload = export.build_site_payload(hall)
    export.write_json(payload, SITE_DATA_DIR / "molecules.json")
    conformers = export.build_conformers_payload(hall)
    export.write_json(conformers, SITE_DATA_DIR / "conformers.json")

    if not catalog_candidates and not mutant_candidates:
        print("\n[run_pipeline] catalogue épuisé et aucun mutant produit "
              "(hall of fame probablement vide — le premier run doit passer par le catalogue).")

    print(f"\n[run_pipeline] {len(passing)} nouveau(x) candidat(s) conforme(s) au TPP "
          f"sur ce run ; hall of fame : {len(hall)} molécule(s) au total.")
    print(f"[run_pipeline] dashboard mis à jour -> {SITE_DATA_DIR}")


if __name__ == "__main__":
    main()
