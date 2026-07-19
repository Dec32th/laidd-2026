from rdkit import Chem
from rdkit.Chem import rdMMPA
from src.tools.replacement_library import get_replacement_candidates

def find_core_and_target(smiles: str, rule_name: str):
    """분자에서 rule_name에 해당하는 문제구조를 담은 조각(target)과
    나머지 뼈대(core)를 찾아서 반환. 못 찾으면 None."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    problem_pattern = Chem.MolFromSmarts(info['problem_smarts'])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    fragments = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)

    for core, chain in fragments:
        if core:
            continue
        parts = chain.split('.')
        if len(parts) != 2:
            continue
        for i, part in enumerate(parts):
            part_mol = Chem.MolFromSmiles(part.replace('[*:1]', '[H]'))
            if part_mol and part_mol.HasSubstructMatch(problem_pattern):
                return {"core": parts[1 - i], "target_removed": part}
    return None


def reassemble_molecule(core_smiles: str, rule_name: str, candidate_idx: int = 0):
    """core의 [*:1] 자리에 replacement_library의 candidate를 붙여 새 분자를 완성."""
    info = get_replacement_candidates(rule_name)
    if info is None or candidate_idx >= len(info['candidates']):
        return None
    candidate = info['candidates'][candidate_idx]

    core_mol = Chem.MolFromSmiles(core_smiles)
    replacement_mol = Chem.MolFromSmiles(f"[*:1]{candidate['smiles']}")
    if core_mol is None or replacement_mol is None:
        return None

    try:
        combined = Chem.molzip(core_mol, replacement_mol)
        new_smiles = Chem.MolToSmiles(combined)
    except Exception:
        return None

    is_valid = Chem.MolFromSmiles(new_smiles) is not None

    return {
        "new_smiles": new_smiles,
        "candidate_used": candidate['name'],
        "rationale": candidate['rationale'],
        "is_valid": is_valid,
    }


def propose_fix(smiles: str, rule_name: str, candidate_idx: int = 0):
    """전체 파이프라인: 문제구조 위치 찾기 -> 치환 후보로 재조립까지 한번에 실행."""
    located = find_core_and_target(smiles, rule_name)
    if located is None:
        return None
    return reassemble_molecule(located['core'], rule_name, candidate_idx)
