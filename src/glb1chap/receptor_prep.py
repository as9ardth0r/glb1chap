"""
glb1chap.receptor_prep
========================
Repérage et extraction du ligand co-cristallisé dans un fichier PDB brut
(téléchargé depuis RCSB) : sert à dériver automatiquement la boîte de
recherche AutoDock Vina autour de sa position réelle plutôt que de deviner
des coordonnées à la main.

Utilisé pour préparer 3THD (bêta-galactosidase humaine, GLB1, en complexe
avec le DGJ) — voir scripts/prepare_receptor.py.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

_IGNORED_RESNAMES = {
    "HOH", "WAT", "DOD",
    "NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "FE2", "CO", "NI",
    "CU", "CU1", "CD", "HG", "BA", "CS", "LI", "RB", "SR", "AL", "PB",
    "AG", "PT", "AU", "GD", "SM", "YB", "LA", "CE", "TB", "IOD", "BR",
    "GOL", "EDO", "PEG", "PG4", "1PE", "2PE", "MPD", "DMS", "SO4",
    "PO4", "ACT", "TRS", "BME", "MRD", "IPA", "FMT", "CIT", "IMD",
    "EPE", "BTB", "BCT", "NO3", "ACY", "CAC", "MES", "BEZ", "PGE",
    "P6G", "1PG", "DIO", "TAM", "UNK", "NAG", "MAN", "BMA", "FUC",
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT",
}


def _read_hetatm_lines(pdb_path: str | Path) -> list[str]:
    lines = Path(pdb_path).read_text(errors="replace").splitlines()
    return [ln for ln in lines if ln.startswith("HETATM")]


def find_cocrystallized_ligand(pdb_path: str | Path) -> tuple[str, str, str, int] | None:
    """Scanne un fichier PDB brut et retourne le ligand co-cristallisé le
    plus probable (chain_id, res_name, res_seq, n_atoms) — dans 3THD, ce
    sera le DGJ. Ignore eau, ions, additifs de cristallisation et glycanes
    de N-glycosylation (fréquents sur une glycoprotéine lysosomale comme
    GLB1)."""
    groups: dict[tuple[str, str, str], int] = defaultdict(int)
    for line in _read_hetatm_lines(pdb_path):
        res_name = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq = line[22:26].strip()
        if res_name in _IGNORED_RESNAMES:
            continue
        groups[(chain_id, res_name, res_seq)] += 1

    if not groups:
        return None

    (chain_id, res_name, res_seq), n_atoms = max(groups.items(), key=lambda kv: kv[1])
    return chain_id, res_name, res_seq, n_atoms


def extract_ligand_pdb(
    pdb_path: str | Path,
    chain_id: str,
    res_name: str,
    res_seq: str,
    output_path: str | Path,
) -> None:
    matching = [
        line for line in _read_hetatm_lines(pdb_path)
        if line[21].strip() == chain_id
        and line[17:20].strip() == res_name
        and line[22:26].strip() == res_seq
    ]
    if not matching:
        raise ValueError(
            f"Aucun atome trouvé pour {res_name} (chaîne {chain_id}, résidu {res_seq}) dans {pdb_path}"
        )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(matching) + "\nEND\n")


def find_all_residue_instances(pdb_path: str | Path, res_name: str) -> list[tuple[str, str]]:
    """3THD est dimérique (2 copies de la protéine dans l'unité asymétrique)
    — le ligand peut apparaître sur chaque copie. Toutes les instances
    doivent être retirées du récepteur préparé."""
    instances: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for line in _read_hetatm_lines(pdb_path):
        if line[17:20].strip() != res_name:
            continue
        key = (line[21].strip(), line[22:26].strip())
        if key not in seen:
            seen.add(key)
            instances.append(key)
    return instances


def format_delete_residues(instances: list[tuple[str, str]]) -> str:
    by_chain: dict[str, list[str]] = {}
    for chain_id, res_seq in instances:
        by_chain.setdefault(chain_id, []).append(res_seq)
    return ",".join(f"{chain}:{','.join(resnums)}" for chain, resnums in by_chain.items())
