"""
glb1chap.generator
====================
Génération combinatoire de candidats chaperons de GLB1 par N-alkylation du
scaffold DGJ (1-désoxygalactonojirimycine) : le noyau piperidine
polyhydroxylé (essentiel à la reconnaissance du site actif, mime le
galactose) est conservé intact, seul le substituant porté par l'azote
endocyclique varie — c'est la stratégie effectivement utilisée pour les
dérivés N-alkylés d'iminosucres (ex. N-butyl-DGJ, miglustat/NB-DNJ) afin de
moduler perméabilité/sélectivité sans toucher au pharmacophore sucre.

Pas de générateur deep-learning : approche rule-based reproductible, sans
dépendance GPU/API externe.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

# Scaffold DGJ avec point d'attache [*:1] sur l'azote endocyclique.
DGJ_SCAFFOLD_SMILES = "OC[C@H]1N([*:1])C[C@H](O)[C@@H](O)[C@H]1O"

# Bibliothèque de substituants N- à explorer (chaînes alkyle de longueurs
# variées, groupes aromatiques/benziliques, chaînes portant un hétéroatome
# pour la solubilité) — inspirée des séries N-alkyl-iminosucres publiées.
R_GROUPS: dict[str, str] = {
    "methyl": "[*:1]C",
    "ethyl": "[*:1]CC",
    "propyl": "[*:1]CCC",
    "butyl": "[*:1]CCCC",
    "pentyl": "[*:1]CCCCC",
    "hexyl": "[*:1]CCCCCC",
    "heptyl": "[*:1]CCCCCCC",
    "octyl": "[*:1]CCCCCCCC",
    "nonyl": "[*:1]CCCCCCCCC",
    "benzyl": "[*:1]Cc1ccccc1",
    "phenethyl": "[*:1]CCc1ccccc1",
    "cyclohexylmethyl": "[*:1]CC1CCCCC1",
    "adamantylmethyl": "[*:1]CC12CC3CC(CC(C3)C1)C2",
    "hydroxyethyl": "[*:1]CCO",
    "carboxypentyl": "[*:1]CCCCCC(=O)O",
    "aminoethyl": "[*:1]CCN",
    "fluorobenzyl": "[*:1]Cc1ccc(F)cc1",
    "methoxybenzyl": "[*:1]Cc1ccc(OC)cc1",
    "morpholinoethyl": "[*:1]CCN1CCOCC1",
    "piperidinylethyl": "[*:1]CCN1CCCCC1",
}


@dataclass
class Recipe:
    scaffold: str
    r_group: str


def _attach(scaffold_smiles: str, r_group_smiles: str) -> str | None:
    """Combine le scaffold (point d'attache [*:1]) et un R-group (même
    label [*:1]) via une réaction RDKit, retourne le SMILES canonique du
    produit, ou None si la combinaison échoue."""
    scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
    r_mol = Chem.MolFromSmiles(r_group_smiles)
    if scaffold_mol is None or r_mol is None:
        return None
    try:
        combined = Chem.molzip(scaffold_mol, r_mol)
    except Exception:
        return None
    if combined is None:
        return None
    try:
        Chem.SanitizeMol(combined)
    except Exception:
        return None
    return Chem.MolToSmiles(combined)


def generate_analog(r_group_name: str) -> tuple[str, Recipe] | None:
    """Génère l'analogue N-substitué correspondant à un nom de R-group de
    `R_GROUPS`. Retourne (smiles, recipe) ou None si la génération échoue."""
    r_group_smiles = R_GROUPS.get(r_group_name)
    if r_group_smiles is None:
        return None
    smiles = _attach(DGJ_SCAFFOLD_SMILES, r_group_smiles)
    if smiles is None:
        return None
    return smiles, Recipe(scaffold="DGJ", r_group=r_group_name)


def generate_all() -> list[tuple[str, Recipe]]:
    """Énumère tous les analogues N-substitués du catalogue R_GROUPS."""
    results = []
    for name in R_GROUPS:
        analog = generate_analog(name)
        if analog is not None:
            results.append(analog)
    return results


def mol_id_for(smiles: str, recipe: Recipe, date: str) -> str:
    digest = hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:8]
    return f"dgj_{recipe.r_group}_{date}_{digest}"
