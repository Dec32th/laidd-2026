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

_DUPLICATE_RULE_MAP = {
    "catechol_A(92)": "catechol",
    "diazo_group": "azo_A(324)",
    "oxime_1": "imine_1_oxime",
    "chinone_1": "quinone_A(370)",
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
    PAINS/BRENK가 동일하거나 부분적으로 겹치는 구조를 서로 다른 이름/원자
    범위로 중복 보고하는 경우(예: catechol_A(92)==catechol,
    diazo_group==azo_A(324), oxime_1이 imine_1_oxime과 원자 하나 차이로
    겹침, chinone_1==quinone_A(370)), 같은 rule_name에 원자 인덱스가
    하나라도 겹치면 중복으로 간주해 제거한다(완전히 동일한 인덱스일
    필요는 없음).
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
                continue
            elif rule_name in _DUPLICATE_RULE_MAP:
                rule_name = _DUPLICATE_RULE_MAP[rule_name]

            is_duplicate = any(
                r['rule_name'] == rule_name and set(r['atom_indices']) & set(atom_indices)
                for r in results
            )
            if is_duplicate:
                continue

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
