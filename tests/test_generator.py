import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem

from glb1chap.generator import R_GROUPS, generate_all, generate_analog, generate_from_fragment, mutate_fragment
from glb1chap.properties import compute_descriptors
from glb1chap.filters import enrich_and_filter, has_basic_amine
from glb1chap import evolve
from glb1chap.hall_of_fame import elite_parents, merge_into_hall_of_fame


def test_all_r_groups_generate_valid_molecules():
    results = generate_all()
    assert len(results) == len(R_GROUPS)
    for smiles, recipe in results:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"SMILES invalide pour {recipe.r_group}: {smiles}"


def test_generated_analogs_preserve_polyol_core():
    smiles, _ = generate_analog("butyl")
    mol = Chem.MolFromSmiles(smiles)
    # 4 groupes hydroxyle attendus (hérités du scaffold DGJ, non touchés par la décoration)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    assert len(mol.GetSubstructMatches(oh_pattern)) == 4


def test_unknown_r_group_returns_none():
    assert generate_analog("does_not_exist") is None


def test_basic_amine_present_on_all_analogs():
    for smiles, recipe in generate_all():
        mol = Chem.MolFromSmiles(smiles)
        assert has_basic_amine(mol), f"amine basique manquante pour {recipe.r_group}"


def test_long_chain_analogs_fail_tpp_on_brenk_alert():
    smiles, recipe = generate_analog("octyl")
    record = compute_descriptors(smiles, "test_octyl")
    record = enrich_and_filter(record)
    assert record.tpp_pass is False
    assert "Brenk" in record.notes


def test_mutation_never_touches_core_polyol():
    import random
    rng = random.Random(7)
    frag = R_GROUPS["butyl"]
    for _ in range(10):
        mutated = mutate_fragment(frag, rng)
        if mutated is not None:
            frag = mutated
    smiles, _ = generate_from_fragment(frag)
    mol = Chem.MolFromSmiles(smiles)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    assert len(mol.GetSubstructMatches(oh_pattern)) == 4
    assert has_basic_amine(mol)


def test_evolve_falls_back_to_mutation_once_catalog_exhausted():
    explored = {smiles for smiles, _ in generate_all()}
    assert evolve.candidates_from_catalog(explored) is not None  # toujours la liste complète, filtrage fait par run_pipeline

    # hall of fame factice avec un parent exploitable
    smiles, recipe = generate_analog("butyl")
    record = compute_descriptors(smiles, "dgj_butyl_test")
    record.recipe = {"scaffold": recipe.scaffold, "r_group": recipe.r_group,
                      "r_group_smiles": recipe.r_group_smiles, "parent_id": None}
    record = enrich_and_filter(record)
    hall = merge_into_hall_of_fame([], [record], "2026-01-01")
    assert elite_parents(hall, n=5)

    mutants = evolve.candidates_from_mutation(hall, explored, n_mutants=5, n_parents=5, seed=1)
    assert len(mutants) > 0
    for smiles, mrecipe in mutants:
        assert smiles not in explored
        assert mrecipe.parent_id == "dgj_butyl_test"
