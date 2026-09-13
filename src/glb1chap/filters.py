"""
glb1chap.filters
=================
Filtres de criblage pour des candidats chaperons de GLB1 (bêta-galactosidase
lysosomale) : profil cible (TPP) adapté aux iminosucres, alertes structurales
PAINS/Brenk (catalogues intégrés RDKit), règle de Lipinski, et un score de
synthétisabilité (SA score) simplifié.

Profil TPP justifié par le mécanisme d'action des chaperons pharmacologiques
de type iminosucre (ex. DGJ/migalastat) : petite molécule polaire, portant
une amine secondaire basique protonable (mime la charge de l'état de
transition de l'hydrolyse enzymatique), suffisamment polaire pour occuper le
site actif mais pas au point de perdre toute perméabilité cellulaire.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import FilterCatalog, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from .properties import MoleculeRecord

# Amine secondaire/tertiaire aliphatique protonable (typique du cycle
# piperidine des iminosucres) — absente => pas de mécanisme chaperon plausible.
_BASIC_AMINE_SMARTS = Chem.MolFromSmarts("[NX3;!$(N-C=[O,N,S]);!$(N-a)]")

def _build_catalog(catalog_enum) -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(catalog_enum)
    return FilterCatalog.FilterCatalog(params)


_PAINS_CATALOG = _build_catalog(FilterCatalogParams.FilterCatalogs.PAINS)
_BRENK_CATALOG = _build_catalog(FilterCatalogParams.FilterCatalogs.BRENK)


@dataclass
class TPPProfile:
    """Profil cible (Target Product Profile) pour un chaperon de GLB1."""

    mw_range: tuple[float, float] = (140.0, 400.0)
    logp_range: tuple[float, float] = (-3.0, 3.0)
    tpsa_range: tuple[float, float] = (40.0, 110.0)
    hbd_range: tuple[int, int] = (2, 6)
    hba_range: tuple[int, int] = (3, 8)
    max_rotatable_bonds: int = 8
    max_lipinski_violations: int = 1
    require_basic_amine: bool = True


DEFAULT_TPP = TPPProfile()


def has_basic_amine(mol: Chem.Mol) -> bool:
    return mol.HasSubstructMatch(_BASIC_AMINE_SMARTS)


def lipinski_violations(record: MoleculeRecord) -> int:
    violations = 0
    if record.mw > 500:
        violations += 1
    if record.logp > 5:
        violations += 1
    if record.hbd > 5:
        violations += 1
    if record.hba > 10:
        violations += 1
    return violations


def pains_alerts(mol: Chem.Mol) -> list[str]:
    return [entry.GetDescription() for entry in _PAINS_CATALOG.GetMatches(mol)]


def brenk_alerts(mol: Chem.Mol) -> list[str]:
    return [entry.GetDescription() for entry in _BRENK_CATALOG.GetMatches(mol)]


def sa_score_heuristic(mol: Chem.Mol) -> float:
    """Score de synthétisabilité simplifié (1 = facile, 10 = très difficile),
    basé sur la complexité structurale (cycles, stéréocentres, hétéroatomes).
    À remplacer par le SA score officiel RDKit Contrib (Ertl & Schuffenhauer
    2009) si `vendor/sascorer.py` est ajouté au projet — non inclus ici."""
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
    n_hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    n_heavy = mol.GetNumHeavyAtoms() or 1
    score = 1.0 + 0.4 * n_rings + 0.5 * n_stereo + 0.15 * n_hetero + 0.02 * n_heavy
    return round(min(score, 10.0), 2)


def enrich_and_filter(record: MoleculeRecord, tpp: TPPProfile = DEFAULT_TPP) -> MoleculeRecord:
    """Calcule les alertes structurales et le verdict TPP pour un record déjà
    doté de descripteurs de base (voir properties.compute_descriptors), et
    retourne le record enrichi (pains_alerts, toxicity_alerts, sa_score,
    lipinski_violations, tpp_pass, notes)."""
    mol = Chem.MolFromSmiles(record.canonical_smiles)
    if mol is None:
        record.tpp_pass = False
        record.notes = "SMILES canonique invalide"
        return record

    record.pains_alerts = pains_alerts(mol)
    record.toxicity_alerts = brenk_alerts(mol)
    record.sa_score = sa_score_heuristic(mol)
    record.lipinski_violations = lipinski_violations(record)

    reasons = []
    if not (tpp.mw_range[0] <= record.mw <= tpp.mw_range[1]):
        reasons.append(f"MW hors plage ({record.mw})")
    if not (tpp.logp_range[0] <= record.logp <= tpp.logp_range[1]):
        reasons.append(f"LogP hors plage ({record.logp})")
    if not (tpp.tpsa_range[0] <= record.tpsa <= tpp.tpsa_range[1]):
        reasons.append(f"TPSA hors plage ({record.tpsa})")
    if not (tpp.hbd_range[0] <= record.hbd <= tpp.hbd_range[1]):
        reasons.append(f"HBD hors plage ({record.hbd})")
    if not (tpp.hba_range[0] <= record.hba <= tpp.hba_range[1]):
        reasons.append(f"HBA hors plage ({record.hba})")
    if record.rotatable_bonds > tpp.max_rotatable_bonds:
        reasons.append("trop de liaisons rotatables")
    if record.lipinski_violations > tpp.max_lipinski_violations:
        reasons.append("trop de violations de Lipinski")
    if tpp.require_basic_amine and not has_basic_amine(mol):
        reasons.append("pas d'amine basique (mécanisme chaperon peu plausible)")
    if record.pains_alerts:
        reasons.append(f"alerte(s) PAINS : {', '.join(record.pains_alerts)}")
    if record.toxicity_alerts:
        reasons.append(f"alerte(s) Brenk : {', '.join(record.toxicity_alerts)}")

    record.tpp_pass = not reasons
    record.notes = "Conforme au TPP" if record.tpp_pass else "; ".join(reasons)
    return record
