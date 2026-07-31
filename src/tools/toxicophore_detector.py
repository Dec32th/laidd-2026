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


def _refine_thiol1(mol, atom_indices):
    """thiol_1은 FilterCatalog 원본이 음이온 황 원자 1개만 매치하는 넓은
    규칙이므로, 그 황이 붙은 탄소의 나머지 결합을 확인해 세분화한다:
    이웃 탄소가 C=S도 가지면 디티오카바메이트, C=O를 가지면
    티오카르복실산염, 둘 다 아니면 일반형(미지원)으로 분류한다."""
    s_idx = atom_indices[0]
    s_atom = mol.GetAtomWithIdx(s_idx)
    for nbr in s_atom.GetNeighbors():
        for nbr2 in nbr.GetNeighbors():
            if nbr2.GetIdx() == s_idx:
                continue
            bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
            if bond2 is None or bond2.GetBondTypeAsDouble() != 2.0:
                continue
            if nbr2.GetSymbol() == 'S':
                return "thiol_1_dithiocarbamate"
            if nbr2.GetSymbol() == 'O':
                return "thiol_1_thiocarboxylate"
    return "thiol_1_general"


def detect_toxicophores(smiles: str) -> list[dict]:
    """
    분자의 SMILES를 받아, FilterCatalog(PAINS+BRENK)에 매치되는
    문제 구조(toxicophore)들을 찾아서 규칙 이름과 해당 원자 인덱스를 반환.
    imine_1은 옥심/구아니딘/일반이민, thiol_1은 디티오카바메이트/
    티오카르복실산염/일반형 하위형으로 세분화하여 반환한다.
    aniline은 FilterCatalog의 단순 [NH2] 탐지 대신, replacement_library의
    확장된 패턴(para-치환 벤젠 포함)을 그대로 사용해 재정의한다.
    PAINS/BRENK가 동일하거나 부분적으로 겹치는 구조를 서로 다른 이름/원자
    범위로 중복 보고하는 경우, 같은 rule_name에 원자 인덱스가 하나라도
    겹치면 중복으로 간주해 제거한다(완전히 동일한 인덱스일 필요는 없음).
    imine_1_general과 isocyanate가 같은 원자(누적이중결합)를 가리키면
    처리 가능한 isocyanate를 우선하고 imine_1_general은 제거한다.
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
            elif rule_name == "thiol_1":
                rule_name = _refine_thiol1(mol, atom_indices)
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

    isocyanate_atom_sets = [set(r['atom_indices']) for r in results if r['rule_name'] == 'isocyanate']
    if isocyanate_atom_sets:
        results = [
            r for r in results
            if not (r['rule_name'] == 'imine_1_general'
                    and any(set(r['atom_indices']) & iso_set for iso_set in isocyanate_atom_sets))
        ]

    return results
