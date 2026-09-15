"""
glb1chap.generator
====================
Génération de candidats chaperons de GLB1 par N-alkylation du scaffold DGJ
(1-désoxygalactonojirimycine) : le noyau piperidine polyhydroxylé (mime le
galactose, essentiel à la reconnaissance du site actif) est TOUJOURS
conservé intact — seul le fragment porté par l'azote endocyclique varie.

Deux mécanismes de génération, combinés par `evolve.py` :
  1. **Exploration par catalogue** (`R_GROUPS`) : un ensemble fini de
     substituants connus/raisonnables, épuisable.
  2. **Mutation atomique** (`mutate_fragment`) : une modification chimique
     locale (ajouter/retirer un halogène, permuter un halogène, ajouter/
     retirer un méthyle) appliquée à un fragment R déjà bon (issu du hall
     of fame). Contrairement à (1), cet espace n'est PAS fini — la mutation
     opère toujours sur le FRAGMENT isolé (jamais sur la molécule complète
     assemblée avec le noyau DGJ), ce qui garantit structurellement que le
     pharmacophore polyol/amine n'est jamais touché, sans avoir besoin de
     le repérer dynamiquement dans la molécule finale.
"""

from __future__ import annotations

import hashlib
import random
from dataclasses import dataclass

from rdkit import Chem

# Scaffold DGJ avec point d'attache [*:1] sur l'azote endocyclique.
DGJ_SCAFFOLD_SMILES = "OC[C@H]1N([*:1])C[C@H](O)[C@@H](O)[C@H]1O"

# Catalogue de substituants N- à explorer : chaînes alkyle (droites et
# ramifiées), cycles alicycliques, groupes (hétéro)aromatiques substitués,
# chaînes fonctionnalisées (OH, NH2, acide, amine cyclique) — inspiré des
# séries N-alkyl-iminosucres publiées (miglustat/NB-DNJ et analogues).
R_GROUPS: dict[str, str] = {
    # chaînes linéaires
    "methyl": "[*:1]C",
    "ethyl": "[*:1]CC",
    "propyl": "[*:1]CCC",
    "butyl": "[*:1]CCCC",
    "pentyl": "[*:1]CCCCC",
    "hexyl": "[*:1]CCCCCC",
    "heptyl": "[*:1]CCCCCCC",
    "octyl": "[*:1]CCCCCCCC",
    "nonyl": "[*:1]CCCCCCCCC",
    # chaînes ramifiées
    "isopropyl": "[*:1]C(C)C",
    "isobutyl": "[*:1]CC(C)C",
    "neopentyl": "[*:1]CC(C)(C)C",
    "2-ethylhexyl": "[*:1]CC(CC)CCCC",
    # cycles alicycliques
    "cyclopropylmethyl": "[*:1]CC1CC1",
    "cyclopentylmethyl": "[*:1]CC1CCCC1",
    "cyclohexylmethyl": "[*:1]CC1CCCCC1",
    "cycloheptylmethyl": "[*:1]CC1CCCCCC1",
    "adamantylmethyl": "[*:1]CC12CC3CC(CC(C3)C1)C2",
    # aromatiques / benzyliques substitués
    "benzyl": "[*:1]Cc1ccccc1",
    "phenethyl": "[*:1]CCc1ccccc1",
    "fluorobenzyl": "[*:1]Cc1ccc(F)cc1",
    "chlorobenzyl": "[*:1]Cc1ccc(Cl)cc1",
    "bromobenzyl": "[*:1]Cc1ccc(Br)cc1",
    "methylbenzyl": "[*:1]Cc1ccc(C)cc1",
    "methoxybenzyl": "[*:1]Cc1ccc(OC)cc1",
    "trifluoromethylbenzyl": "[*:1]Cc1ccc(C(F)(F)F)cc1",
    "naphthylmethyl": "[*:1]Cc1ccc2ccccc2c1",
    "furylmethyl": "[*:1]Cc1ccco1",
    "thienylmethyl": "[*:1]Cc1cccs1",
    "pyridylmethyl": "[*:1]Cc1ccccn1",
    # chaînes fonctionnalisées
    "hydroxyethyl": "[*:1]CCO",
    "hydroxypropyl": "[*:1]CCCO",
    "carboxypentyl": "[*:1]CCCCCC(=O)O",
    "aminoethyl": "[*:1]CCN",
    "cyanoethyl": "[*:1]CCC#N",
    "fluoropropyl": "[*:1]CCCF",
    "difluoroethyl": "[*:1]CC(F)F",
    # amines cycliques (solubilité/biodisponibilité)
    "morpholinoethyl": "[*:1]CCN1CCOCC1",
    "piperidinylethyl": "[*:1]CCN1CCCCC1",
    "pyrrolidinylethyl": "[*:1]CCN1CCCC1",
    "piperazinylethyl": "[*:1]CCN1CCNCC1",
}

# Halogènes utilisés par la mutation (add_halogen / swap_halogen).
_HALOGENS = ("F", "Cl", "Br")
_HALOGEN_ATOMIC_NUMS = {"F": 9, "Cl": 17, "Br": 35}


@dataclass
class Recipe:
    scaffold: str
    r_group: str  # nom du catalogue, ou "mutant" pour un fragment muté
    r_group_smiles: str  # fragment [*:1]... exact utilisé (reproductibilité + mutation future)
    parent_id: str | None = None  # id du parent muté, None sinon


def _attach(scaffold_smiles: str, r_group_smiles: str) -> str | None:
    """Combine le scaffold (point d'attache [*:1]) et un fragment R (même
    label [*:1]) via une réaction RDKit, retourne le SMILES canonique du
    produit, ou None si la combinaison échoue."""
    scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
    r_mol = Chem.MolFromSmiles(r_group_smiles)
    if scaffold_mol is None or r_mol is None:
        return None
    try:
        combined = Chem.molzip(scaffold_mol, r_mol)
        Chem.SanitizeMol(combined)
    except Exception:
        return None
    return Chem.MolToSmiles(combined)


def generate_analog(r_group_name: str) -> tuple[str, Recipe] | None:
    """Génère l'analogue N-substitué correspondant à un nom de R-group du
    catalogue `R_GROUPS`. Retourne (smiles, recipe) ou None si échec."""
    r_group_smiles = R_GROUPS.get(r_group_name)
    if r_group_smiles is None:
        return None
    smiles = _attach(DGJ_SCAFFOLD_SMILES, r_group_smiles)
    if smiles is None:
        return None
    return smiles, Recipe(scaffold="DGJ", r_group=r_group_name, r_group_smiles=r_group_smiles)


def generate_from_fragment(r_group_smiles: str, name: str = "mutant", parent_id: str | None = None) -> tuple[str, Recipe] | None:
    """Génère l'analogue correspondant à un fragment R arbitraire (utilisé
    pour les mutants produits par `mutate_fragment`, absents du catalogue)."""
    smiles = _attach(DGJ_SCAFFOLD_SMILES, r_group_smiles)
    if smiles is None:
        return None
    return smiles, Recipe(scaffold="DGJ", r_group=name, r_group_smiles=r_group_smiles, parent_id=parent_id)


def generate_all() -> list[tuple[str, Recipe]]:
    """Énumère tous les analogues N-substitués du catalogue R_GROUPS."""
    results = []
    for name in R_GROUPS:
        analog = generate_analog(name)
        if analog is not None:
            results.append(analog)
    return results


def mutate_fragment(r_group_smiles: str, rng: random.Random | None = None) -> str | None:
    """Applique UNE modification chimique locale au fragment R (jamais au
    noyau DGJ, puisque le noyau n'existe même pas dans le mol traité ici) :
    ajouter un halogène, en retirer un, en permuter un, ajouter un méthyle,
    ou en retirer un (méthyle terminal). Retourne le nouveau fragment
    SMILES (avec point d'attache [*:1] préservé), ou None si aucune
    mutation n'est applicable ou si le résultat est invalide."""
    rng = rng or random.Random()
    mol = Chem.MolFromSmiles(r_group_smiles)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)

    carbons_with_h = [
        a.GetIdx() for a in rw.GetAtoms()
        if a.GetSymbol() == "C" and a.GetTotalNumHs() > 0
    ]
    halogens = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() in _HALOGENS]
    terminal_methyls = [
        a.GetIdx() for a in rw.GetAtoms()
        if a.GetSymbol() == "C" and a.GetDegree() == 1 and a.GetTotalNumHs() == 3
    ]

    moves = []
    if carbons_with_h:
        moves += ["add_halogen"] * 2 + ["add_methyl"] * 2
    if halogens:
        moves += ["remove_halogen", "swap_halogen"]
    if terminal_methyls:
        moves += ["remove_methyl"]

    if not moves:
        return None
    move = rng.choice(moves)

    try:
        if move == "add_halogen":
            idx = rng.choice(carbons_with_h)
            new_idx = rw.AddAtom(Chem.Atom(rng.choice(_HALOGENS)))
            rw.AddBond(idx, new_idx, Chem.BondType.SINGLE)
        elif move == "add_methyl":
            idx = rng.choice(carbons_with_h)
            new_idx = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(idx, new_idx, Chem.BondType.SINGLE)
        elif move == "remove_halogen":
            idx = rng.choice(halogens)
            rw.RemoveAtom(idx)
        elif move == "swap_halogen":
            idx = rng.choice(halogens)
            current = rw.GetAtomWithIdx(idx).GetSymbol()
            choices = [h for h in _HALOGENS if h != current]
            rw.GetAtomWithIdx(idx).SetAtomicNum(_HALOGEN_ATOMIC_NUMS[rng.choice(choices)])
        elif move == "remove_methyl":
            idx = rng.choice(terminal_methyls)
            rw.RemoveAtom(idx)

        mutated = rw.GetMol()
        Chem.SanitizeMol(mutated)
        return Chem.MolToSmiles(mutated)
    except Exception:
        return None


def mol_id_for(smiles: str, recipe: Recipe, date: str) -> str:
    digest = hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:8]
    return f"dgj_{recipe.r_group}_{date}_{digest}"
