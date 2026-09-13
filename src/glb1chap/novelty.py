"""
glb1chap.novelty
===============
Vérification de nouveauté : la molécule existe-t-elle déjà dans une base
publique (PubChem) ? Comparaison EXACTE par InChIKey.

Étant donné le catalogue de R-groups volontairement restreint (~20 analogues
N-substitués du DGJ), il est très possible — et attendu — que plusieurs
d'entre eux (ex. N-butyl-DGJ) soient déjà des composés publiés/connus :
`is_novel=False` n'est pas un échec, c'est une information utile (littérature
existante à consulter en priorité).
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request

from rdkit import Chem

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def inchikey_for(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


def check_pubchem(inchikey: str, timeout: float = 10.0) -> dict | None:
    """Interroge PubChem PUG-REST par InChIKey exact. Retourne un dict
    {cid, iupac_name} si trouvé, None sinon (y compris en cas d'erreur
    réseau — ne doit jamais interrompre le pipeline)."""
    url = (
        f"{PUBCHEM_BASE}/compound/inchikey/{urllib.parse.quote(inchikey)}"
        f"/property/IUPACName/JSON"
    )
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        return {
            "cid": props[0].get("CID"),
            "iupac_name": props[0].get("IUPACName"),
        }
    except Exception:
        return None
