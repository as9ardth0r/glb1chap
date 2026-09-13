"""
glb1chap.properties
====================
Chargement de SMILES et calcul des descripteurs physico-chimiques
(MW, LogP, TPSA, HBD, HBA, rotatable bonds, QED) via RDKit.

Module générique, indépendant de la cible : toute la suite (filtres,
génération, export) consomme des `MoleculeRecord`.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


@dataclass
class MoleculeRecord:
    """Représentation d'une molécule et de ses descripteurs calculés."""

    id: str
    smiles: str
    canonical_smiles: str
    mw: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    num_rings: int
    num_aromatic_rings: int
    qed: float
    formula: str
    sa_score: Optional[float] = None
    lipinski_violations: Optional[int] = None
    pains_alerts: list[str] = field(default_factory=list)
    toxicity_alerts: list[str] = field(default_factory=list)  # alertes BRENK
    tpp_pass: Optional[bool] = None
    notes: str = ""
    source: str = "generated"
    recipe: Optional[dict] = None  # {"scaffold": "DGJ", "n_substituent": ...}
    first_seen: Optional[str] = None
    fitness: Optional[float] = None
    docking_score: Optional[float] = None  # kcal/mol (AutoDock Vina) contre GLB1 (3THD)
    is_novel: Optional[bool] = None
    pubchem_cid: Optional[int] = None
    chembl_id: Optional[str] = None
    chemical_name: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


class InvalidSMILESError(ValueError):
    """Levée quand une chaîne SMILES ne peut pas être parsée par RDKit."""


def load_smiles(source: str | Path) -> list[tuple[str, str]]:
    """
    Charge une liste de (id, smiles) depuis :
      - un fichier .smi / .txt  (une entrée par ligne : "SMILES id [tag]")
      - un fichier .csv         (colonnes 'id' et 'smiles')
      - une simple chaîne SMILES unique
    """
    path = Path(source) if isinstance(source, (str, Path)) else None
    if path is not None and path.exists():
        if path.suffix == ".csv":
            with path.open(encoding="utf-8") as f:
                reader = csv.DictReader(f)
                return [(row["id"], row["smiles"]) for row in reader]
        entries = []
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            mol_id = parts[1] if len(parts) > 1 else smiles
            entries.append((mol_id, smiles))
        return entries
    # chaîne SMILES unique
    return [(str(source), str(source))]


def compute_descriptors(smiles: str, mol_id: str) -> MoleculeRecord:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"SMILES invalide : {smiles!r}")

    canonical = Chem.MolToSmiles(mol)
    ri = mol.GetRingInfo()
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    return MoleculeRecord(
        id=mol_id,
        smiles=smiles,
        canonical_smiles=canonical,
        mw=round(Descriptors.MolWt(mol), 2),
        logp=round(Crippen.MolLogP(mol), 2),
        tpsa=round(rdMolDescriptors.CalcTPSA(mol), 2),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        num_rings=ri.NumRings(),
        num_aromatic_rings=num_aromatic_rings,
        qed=round(QED.qed(mol), 3),
        formula=rdMolDescriptors.CalcMolFormula(mol),
    )


def compute_batch(entries: list[tuple[str, str]]) -> list[MoleculeRecord]:
    """Calcule les descripteurs pour un lot de (id, smiles), en ignorant
    silencieusement les SMILES invalides (loggés, pas levés) pour ne pas
    interrompre un run batch à cause d'une seule molécule mal formée."""
    records = []
    for mol_id, smiles in entries:
        try:
            records.append(compute_descriptors(smiles, mol_id))
        except InvalidSMILESError:
            logger.warning("SMILES invalide ignoré : id=%s smiles=%s", mol_id, smiles)
    return records


def load_records(path: str | Path) -> list[MoleculeRecord]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        raw = json.load(f)
    return [MoleculeRecord(**entry) for entry in raw]


def save_records(records: list[MoleculeRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump([r.to_dict() for r in records], f, indent=2, ensure_ascii=False)
