"""
glb1chap.hall_of_fame
====================
Persistance de l'état entre deux exécutions du pipeline.

Deux fichiers dans data/ :
  - hall_of_fame.json : les meilleures molécules conformes au TPP à ce jour,
                         plafonnées à HALL_OF_FAME_MAX, classées par fitness.
  - explored.json      : SMILES canoniques déjà testés, pour ne jamais
                         retester deux fois le même analogue.
"""

from __future__ import annotations

import json
from pathlib import Path

from .properties import MoleculeRecord

HALL_OF_FAME_MAX = 100  # catalogue de R-groups volontairement petit (~20) :
                          # pas besoin d'un plafond aussi large que pour un
                          # espace combinatoire ouvert.


def fitness(record: MoleculeRecord) -> float:
    """Score composite simple : favorise le docking (si disponible), sinon
    retombe sur le QED — le docking étant l'évaluation la plus pertinente
    de l'affinité pour le site actif de GLB1."""
    if record.docking_score is not None:
        # docking_score est négatif (kcal/mol) ; on veut un score croissant
        # avec l'affinité, normalisé grossièrement sur une plage -12..0.
        return round(max(0.0, min(1.0, (-record.docking_score) / 12.0)), 4)
    return record.qed


def load_hall_of_fame(path: str | Path) -> list[MoleculeRecord]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        raw = json.load(f)
    return [MoleculeRecord(**entry) for entry in raw]


def save_hall_of_fame(records: list[MoleculeRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump([r.to_dict() for r in records], f, indent=2, ensure_ascii=False)


def load_explored(path: str | Path) -> set[str]:
    path = Path(path)
    if not path.exists():
        return set()
    with path.open(encoding="utf-8") as f:
        return set(json.load(f))


def save_explored(explored: set[str], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(sorted(explored), f, indent=2, ensure_ascii=False)


def merge_into_hall_of_fame(
    hall: list[MoleculeRecord],
    new_passing: list[MoleculeRecord],
    today: str,
    max_size: int = HALL_OF_FAME_MAX,
) -> list[MoleculeRecord]:
    by_canonical: dict[str, MoleculeRecord] = {r.canonical_smiles: r for r in hall}
    for r in new_passing:
        r.fitness = fitness(r)
        if r.canonical_smiles not in by_canonical:
            r.first_seen = today
            by_canonical[r.canonical_smiles] = r

    merged = list(by_canonical.values())
    for r in merged:
        if r.fitness is None:
            r.fitness = fitness(r)
    merged.sort(key=lambda r: r.fitness, reverse=True)
    return merged[:max_size]
