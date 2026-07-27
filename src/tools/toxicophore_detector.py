from rdkit import Chem
from rdkit.Chem import FilterCatalog

def _build_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)

_catalog = _build_catalog()
_oxime_pattern = Chem.MolFromSmarts("C=N[OX2H1]")
_para_aniline_pattern = Chem.MolFromSmarts("[NH2]c1ccc([#6,#7,#8,#16])cc1")


def _refine_imine1(mol, atom_indices):
    """imine_1은 옥심(C=N-OH)과 일반 이민(C=N-R)을 모두 포함하는 넓은 카테고리이므로,
    실제 매치된 부분이 옥심 패턴을 포함하는지 확인해 이름을 세분화한다."""
    if mol.HasSubstructMatch(_oxime_pattern):
        matches = mol.GetSubstructMatches(_oxime_pattern)
        for match in matches:
            if set(match) & set(atom_indices):
                return "imine_1_oxime"
    return "imine_1_general"


def detect_toxicophores(smiles: str) -> list[dict]:
    """
    분자의 SMILES를 받아, FilterCatalog(PAINS+BRENK)에 매치되는
    문제 구조(toxicophore)들을 찾아서 규칙 이름과 해당 원자 인덱스를 반환.
    imine_1은 옥심/일반이민 하위형으로 세분화하여 반환한다.
    추가로, para-치환 아닐린(문헌 기반 커스텀 규칙, FilterCatalog 항목 아님)을
    별도로 탐지하여 aniline과 함께 반환한다.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    results = []
    for entry in _catalog.GetMatches(mol):
        for fm in entry.GetFilterMatches(mol):
            atom_indices = sorted(set(mol_idx for _, mol_idx in fm.atomPairs))
            rule_name = entry.GetDescription()

            if rule_name == "imine_1":
                rule_name = _refine_imine1(mol, atom_indices)

            results.append({
                "rule_name": rule_name,
                "atom_indices": atom_indices,
            })

    # 커스텀 규칙: para-치환 아닐린 (FilterCatalog에 없는, 문헌 기반 자체 추가 규칙)
    if mol.HasSubstructMatch(_para_aniline_pattern):
        matches = mol.GetSubstructMatches(_para_aniline_pattern)
        for match in matches:
            atom_indices = sorted(set(match))
            results.append({
                "rule_name": "aniline_ring_bcp",
                "atom_indices": atom_indices,
            })

    return results
