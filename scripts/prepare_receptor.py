#!/usr/bin/env python3
"""Prépare la structure réceptrice GLB1 pour le docking (docking.py) :
télécharge 3THD depuis RCSB (bêta-galactosidase humaine + DGJ co-cristallisé),
repère le DGJ pour dériver automatiquement la boîte de recherche autour de
sa position réelle, puis appelle `mk_prepare_receptor.py` (meeko) pour
produire le PDBQT et la configuration de boîte Vina.

Nécessite un accès réseau vers files.rcsb.org.

Usage :
    python scripts/prepare_receptor.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from glb1chap.receptor_prep import (
    extract_ligand_pdb,
    find_cocrystallized_ligand,
    find_all_residue_instances,
    format_delete_residues,
)  # noqa: E402

TARGET_NAME = "GLB1"
PDB_ID = "3THD"  # bêta-galactosidase humaine + DGJ (1.79 Å, ligand co-cristallisé)

REPO_ROOT = Path(__file__).resolve().parent.parent
RECEPTOR_DIR = REPO_ROOT / "data" / "receptor" / TARGET_NAME


def fetch_pdb(pdb_id: str, out_path: Path) -> None:
    import requests
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    out_path.write_text(response.text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--padding", type=float, default=4.0,
                         help="marge (Å) ajoutée autour du DGJ pour la boîte de recherche")
    parser.add_argument("--allow-bad-residues", action=argparse.BooleanOptionalAction, default=True,
                         help="supprime les résidus sans template plutôt que de stopper sur erreur")
    parser.add_argument("--default-altloc", type=str, default="A")
    args = parser.parse_args()

    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    raw_path = RECEPTOR_DIR / f"{PDB_ID}.pdb"
    print(f"[prepare_receptor] téléchargement {PDB_ID} -> {raw_path}")
    fetch_pdb(PDB_ID, raw_path)

    ligand_hit = find_cocrystallized_ligand(raw_path)
    output_basename = str(RECEPTOR_DIR / TARGET_NAME.lower())
    cmd = [
        "mk_prepare_receptor.py",
        "--read_pdb", str(raw_path),
        "-o", output_basename,
        "-p", "-v",
        "--padding", str(args.padding),
        "--default_altloc", args.default_altloc,
    ]
    if args.allow_bad_residues:
        cmd.append("-a")

    if ligand_hit is None:
        print("[prepare_receptor] ATTENTION — aucun ligand co-cristallisé trouvé : "
              "boîte à définir manuellement (--box_center / --box_size).")
    else:
        chain_name, res_name, res_seq, n_atoms = ligand_hit
        print(f"[prepare_receptor] ligand co-cristallisé repéré : "
              f"{res_name} (chaîne {chain_name}, résidu {res_seq}, {n_atoms} atomes) — attendu : DGJ")
        ligand_path = RECEPTOR_DIR / f"{TARGET_NAME.lower()}_ref_ligand.pdb"
        extract_ligand_pdb(raw_path, chain_name, res_name, res_seq, ligand_path)
        cmd += ["--box_enveloping", str(ligand_path)]

        instances = find_all_residue_instances(raw_path, res_name)
        delete_arg = format_delete_residues(instances)
        cmd += ["--delete_residues", delete_arg]
        print(f"[prepare_receptor] exclusion du récepteur : {res_name} "
              f"({len(instances)} copie(s), structure dimérique : {delete_arg})")

    print(f"[prepare_receptor] exécution : {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[prepare_receptor] {TARGET_NAME} prêt : "
          f"{output_basename}.pdbqt + {output_basename}_box.txt")


if __name__ == "__main__":
    main()
