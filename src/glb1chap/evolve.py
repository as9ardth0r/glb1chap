"""
glb1chap.evolve
==============
Combine les deux mécanismes de génération pour produire les candidats d'un
run :
  1. **Catalogue** : tous les analogues de `generator.R_GROUPS` pas encore
     dans `explored` (fini — s'épuise après le premier run complet).
  2. **Mutation** : à partir des `n_parents` meilleures molécules du hall
     of fame, applique une mutation atomique aléatoire à leur fragment R
     (`generator.mutate_fragment`), jusqu'à produire `n_mutants` candidats
     inédits. Espace ouvert — c'est ce mécanisme qui permet au pipeline de
     continuer à générer indéfiniment une fois le catalogue épuisé.

Le noyau DGJ n'est jamais manipulé directement : la mutation opère sur le
fragment R isolé, avant assemblage (`generator.generate_from_fragment`).
"""

from __future__ import annotations

import random

from .generator import (
    R_GROUPS,
    generate_analog,
    generate_from_fragment,
    mutate_fragment,
)
from .hall_of_fame import elite_parents
from .properties import MoleculeRecord


def candidates_from_catalog(explored: set[str]) -> list[tuple[str, "Recipe"]]:
    """Analogues du catalogue dont le SMILES canonique n'est pas déjà dans
    `explored`. Le filtrage par canonique se fait ici sur le SMILES brut
    (pas encore canonicalisé) — `run_pipeline.py` refait le test après
    calcul des descripteurs, qui est la source de vérité."""
    results = []
    for name in R_GROUPS:
        analog = generate_analog(name)
        if analog is not None:
            results.append(analog)
    return results


def candidates_from_mutation(
    hall: list[MoleculeRecord],
    explored: set[str],
    n_mutants: int = 20,
    n_parents: int = 5,
    max_attempts: int = 200,
    seed: int | None = None,
) -> list[tuple[str, "Recipe"]]:
    """Génère jusqu'à `n_mutants` candidats inédits par mutation atomique
    des meilleurs parents du hall of fame. Retourne une liste possiblement
    plus courte que `n_mutants` si le hall of fame est vide (premier run,
    avant qu'aucune molécule n'ait encore été retenue) ou si les tentatives
    de mutation s'épuisent sans produire de nouveauté."""
    rng = random.Random(seed)
    parents = elite_parents(hall, n=n_parents)
    if not parents:
        return []

    results = []
    seen_this_run: set[str] = set()
    attempts = 0
    while len(results) < n_mutants and attempts < max_attempts:
        attempts += 1
        parent = rng.choice(parents)
        parent_fragment = parent.recipe["r_group_smiles"]
        mutated_fragment = mutate_fragment(parent_fragment, rng)
        if mutated_fragment is None:
            continue
        analog = generate_from_fragment(mutated_fragment, name="mutant", parent_id=parent.id)
        if analog is None:
            continue
        smiles, recipe = analog
        if smiles in explored or smiles in seen_this_run:
            continue
        seen_this_run.add(smiles)
        results.append(analog)
    return results
