from rdkit import Chem
from rdkit.Chem import rdMMPA
from src.tools.replacement_library import get_replacement_candidates


def find_core_and_target(smiles: str, rule_name: str):
    """분자에서 rule_name에 해당하는 문제구조를 담은 조각(target)과
    나머지 뼈대(core)를 찾아서 반환. 패턴 크기와 정확히 일치하는 조각만 인정.
    못 찾으면 None (차선책으로 얼버무리지 않음)."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    problem_pattern = Chem.MolFromSmarts(info['problem_smarts'])
    pattern_size = problem_pattern.GetNumAtoms()

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
            if part_mol is None or not part_mol.HasSubstructMatch(problem_pattern):
                continue
            frag_heavy_atoms = part_mol.GetNumHeavyAtoms()
            if frag_heavy_atoms == pattern_size:
                return {"core": parts[1 - i], "target_removed": part}
    return None


def reassemble_molecule(core_smiles: str, rule_name: str, candidate_idx: int = 0):
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
    located = find_core_and_target(smiles, rule_name)
    if located is None:
        return None
    return reassemble_molecule(located['core'], rule_name, candidate_idx)


def canonicalize(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def iterative_fix_loop(smiles: str, max_iterations: int = 10, candidate_idx: int = 0,
                        llm_client=None, llm_model=None):
    """진단->치환->재평가를 반복.
    llm_client가 주어지면: 어떤 문제부터 고칠지 + 어떤 후보를 쓸지 둘 다 LLM이 판단.
    없으면: 리스트 순서(known_problems[0]) + candidate_idx 고정값 사용."""
    from src.tools.toxicophore_detector import detect_toxicophores
    from src.tools.agent import ask_llm_which_problem_to_fix, ask_llm_which_candidate_to_use

    current = canonicalize(smiles)
    seen = {current}
    history = [{"step": 0, "smiles": current}]
    skipped_rules = []

    for step in range(1, max_iterations + 1):
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

        if llm_client is not None:
            problem_decision = ask_llm_which_problem_to_fix(llm_client, llm_model, current, problems)
            target_rule = problem_decision['rule_name']
            problem_reason = problem_decision.get('reason', '')

            candidate_decision = ask_llm_which_candidate_to_use(llm_client, llm_model, current, target_rule)
            chosen_candidate_idx = candidate_decision['candidate_idx']
            candidate_reason = candidate_decision.get('reason', '')
        else:
            target_rule = known_problems[0]['rule_name']
            problem_reason = "규칙 기반(리스트 순서대로)"
            chosen_candidate_idx = candidate_idx
            candidate_reason = "규칙 기반(고정 인덱스)"

        fixed = propose_fix(current, target_rule, chosen_candidate_idx)

        if fixed is None or not fixed['is_valid']:
            return {"status": "stuck", "reason": f"'{target_rule}' 치환 실패 (core/target 매칭 실패 또는 재조립 실패)", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        new_current = canonicalize(fixed['new_smiles'])

        if new_current in seen:
            return {"status": "cycle_detected", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}

        seen.add(new_current)
        current = new_current
        history.append({
            "step": step,
            "smiles": current,
            "fixed_rule": target_rule,
            "problem_reason": problem_reason,
            "candidate_used": fixed['candidate_used'],
            "candidate_reason": candidate_reason,
        })

    return {"status": "max_iterations_reached", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}
