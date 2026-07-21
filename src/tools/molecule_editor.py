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
    pattern_size = problem_pattern.GetNumAtoms()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    fragments = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)

    best_match = None
    for core, chain in fragments:
        if core:
            continue
        parts = chain.split('.')
        if len(parts) != 2:
            continue
        for i, part in enumerate(parts):
            part_mol = Chem.MolFromSmiles(part.replace('[*:1]', '[H]'))
            if part_mol is None or not part_mol.HasSubstructMatch(problem_pattern):
                continue
            frag_heavy_atoms = part_mol.GetNumHeavyAtoms()
            if frag_heavy_atoms == pattern_size:
                return {"core": parts[1 - i], "target_removed": part}
            if best_match is None:
                best_match = {"core": parts[1 - i], "target_removed": part}
    return best_match


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


def canonicalize(smiles: str):
    """SMILES를 canonical(정규) 형태로 변환. 파싱 실패 시 None."""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def iterative_fix_loop(smiles: str, max_iterations: int = 10, candidate_idx: int = 0):
    """진단->치환->재평가를 반복. 성공/실패/순환/미지의 규칙 등으로 종료."""
    current = canonicalize(smiles)
    seen = {current}
    history = [{"step": 0, "smiles": current}]
    skipped_rules = []

    for step in range(1, max_iterations + 1):
        # 순환참조 방지를 위해 여기서 import (같은 파일 내 함수는 아래에서 직접 씀)
        from src.tools.toxicophore_detector import detect_toxicophores
        problems = detect_toxicophores(current)
        history[-1]["problems"] = problems

        if not problems:
            return {"status": "success", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        known_problems = [p for p in problems if get_replacement_candidates(p['rule_name']) is not None]
        unknown_problems = [p for p in problems if p not in known_problems]

        for p in unknown_problems:
            if p['rule_name'] not in skipped_rules:
                skipped_rules.append(p['rule_name'])

        if not known_problems:
            return {"status": "no_known_fix", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        target_rule = known_problems[0]['rule_name']
        fixed = propose_fix(current, target_rule, candidate_idx)

        if fixed is None or not fixed['is_valid']:
            return {"status": "stuck", "reason": f"'{target_rule}' 치환 실패", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        new_current = canonicalize(fixed['new_smiles'])

        if new_current in seen:
            return {"status": "cycle_detected", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        seen.add(new_current)
        current = new_current
        history.append({
            "step": step,
            "smiles": current,
            "fixed_rule": target_rule,
            "candidate_used": fixed['candidate_used'],
        })

    return {"status": "max_iterations_reached", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}
