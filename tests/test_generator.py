import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem

from glb1chap.generator import R_GROUPS, generate_all, generate_analog
from glb1chap.properties import compute_descriptors
from glb1chap.filters import enrich_and_filter, has_basic_amine


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
