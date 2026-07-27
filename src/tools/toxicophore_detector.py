from rdkit import Chem
from rdkit.Chem import FilterCatalog

def _build_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)

_catalog = _build_catalog()
_oxime_pattern = Chem.MolFromSmarts("C=N[OX2H1]")


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
    aniline은 FilterCatalog의 단순 [NH2] 탐지 대신, replacement_library의
    확장된 패턴(para-치환 벤젠 포함)을 그대로 사용해 재정의한다.
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
            elif rule_name == "aniline":
                continue  # 아래에서 확장된 패턴으로 다시 탐지하므로 원본은 건너뜀

            results.append({
                "rule_name": rule_name,
                "atom_indices": atom_indices,
            })

    # aniline: 원래 FilterCatalog의 단순 [NH2] 대신, replacement_library와
    # 정확히 일치하는 확장된 SMARTS로 재탐지 (아민만 vs 고리 전체, 두 candidate가
    # 모두 이 하나의 매치 위에서 작동하도록 통일)
    from src.tools.replacement_library import get_replacement_candidates
    aniline_info = get_replacement_candidates("aniline")
    if aniline_info:
        aniline_pattern = Chem.MolFromSmarts(aniline_info["problem_smarts"])
        if mol.HasSubstructMatch(aniline_pattern):
            matches = mol.GetSubstructMatches(aniline_pattern)
            for match in matches:
                atom_indices = sorted(set(match))
                results.append({
                    "rule_name": "aniline",
                    "atom_indices": atom_indices,
                })

    return results
