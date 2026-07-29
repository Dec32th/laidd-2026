from rdkit import Chem
from rdkit.Chem import FilterCatalog

def _build_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)

_catalog = _build_catalog()
_oxime_pattern = Chem.MolFromSmarts("C=N[OX2H1]")
_guanidine_pattern = Chem.MolFromSmarts("[$(C(N)(N)=N)]")

# PAINS와 BRENK 양쪽에 동일 화학구조를 잡는 중복 규칙명이 있는 경우,
# 우리 라이브러리 기준 이름으로 통일 (동일 원자 인덱스로 확인된 것만)
_DUPLICATE_RULE_MAP = {
    "catechol_A(92)": "catechol",
}


def _refine_imine1(mol, atom_indices):
    """imine_1은 옥심(C=N-OH), 구아니딘(N-C(=N)-N), 일반 이민(C=N-R)을
    모두 포함하는 넓은 카테고리이므로, 실제 매치 부분의 화학적 맥락을
    확인해 이름을 세분화한다."""
    if mol.HasSubstructMatch(_oxime_pattern):
        matches = mol.GetSubstructMatches(_oxime_pattern)
        for match in matches:
            if set(match) & set(atom_indices):
                return "imine_1_oxime"
    if mol.HasSubstructMatch(_guanidine_pattern):
        matches = mol.GetSubstructMatches(_guanidine_pattern)
        for match in matches:
            if set(match) & set(atom_indices):
                return "imine_1_guanidine"
    return "imine_1_general"


def detect_toxicophores(smiles: str) -> list[dict]:
    """
    분자의 SMILES를 받아, FilterCatalog(PAINS+BRENK)에 매치되는
    문제 구조(toxicophore)들을 찾아서 규칙 이름과 해당 원자 인덱스를 반환.
    imine_1은 옥심/구아니딘/일반이민 하위형으로 세분화하여 반환한다.
    aniline은 FilterCatalog의 단순 [NH2] 탐지 대신, replacement_library의
    확장된 패턴(para-치환 벤젠 포함)을 그대로 사용해 재정의한다.
    PAINS/BRENK가 동일 원자를 서로 다른 이름으로 중복 보고하는 경우
    (예: catechol_A(92) == catechol), 라이브러리 기준 이름으로 통일하고
    중복 항목은 제거한다.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    results = []
    seen_entries = set()

    for entry in _catalog.GetMatches(mol):
        for fm in entry.GetFilterMatches(mol):
            atom_indices = sorted(set(mol_idx for _, mol_idx in fm.atomPairs))
            rule_name = entry.GetDescription()

            if rule_name == "imine_1":
                rule_name = _refine_imine1(mol, atom_indices)
            elif rule_name == "aniline":
                continue
            elif rule_name in _DUPLICATE_RULE_MAP:
                rule_name = _DUPLICATE_RULE_MAP[rule_name]

            dedup_key = (rule_name, tuple(atom_indices))
            if dedup_key in seen_entries:
                continue
            seen_entries.add(dedup_key)

            results.append({
                "rule_name": rule_name,
                "atom_indices": atom_indices,
            })

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
