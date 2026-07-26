from rdkit import Chem
from rdkit.Chem import rdMMPA
from src.tools.replacement_library import get_replacement_candidates


def _check_and_match(part_smiles, problem_pattern, pattern_size):
    """조각이 problem_pattern과 정확한 크기로 매치되는지 확인."""
    part_mol = Chem.MolFromSmiles(part_smiles.replace('[*:1]', 'C').replace('[*:2]', 'C'))
    if part_mol is None or not part_mol.HasSubstructMatch(problem_pattern):
        return False
    n_attachment = part_smiles.count('[*:')
    return part_mol.GetNumHeavyAtoms() - n_attachment == pattern_size


def find_core_and_target(smiles: str, rule_name: str):
    """분자에서 rule_name에 해당하는 문제구조를 담은 조각(target)과
    나머지 뼈대(core)를 찾아서 반환.
    1단계(maxCuts=1)로 단순 분리를 먼저 시도하고,
    실패하면 2단계(maxCuts=2)로 고리 인접 작용기 분리를 시도한다."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    problem_pattern = Chem.MolFromSmarts(info['problem_smarts'])
    pattern_size = problem_pattern.GetNumAtoms()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # --- Case A: 단순 2조각 분리 (maxCuts=1) ---
    fragments1 = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)
    for core, chain in fragments1:
        if core:
            continue
        parts = chain.split('.')
        if len(parts) != 2:
            continue
        for i, part in enumerate(parts):
            if _check_and_match(part, problem_pattern, pattern_size):
                return {"core": parts[1 - i], "target_removed": part}

    # --- Case B: 고리 인접 등, core가 남는 2-cut 분리 ---
    fragments2 = rdMMPA.FragmentMol(mol, maxCuts=2, resultsAsMols=False)
    for core, chain in fragments2:
        if not core:
            continue
        chain_parts = chain.split('.')
        if len(chain_parts) != 2:
            continue
        for i, part in enumerate(chain_parts):
            if not _check_and_match(part, problem_pattern, pattern_size):
                continue
            other_chain_part = chain_parts[1 - i]

            target_ap = '[*:1]' if '[*:1]' in part else ('[*:2]' if '[*:2]' in part else None)
            if target_ap is None:
                continue

            core_mol = Chem.MolFromSmiles(core)
            other_mol = Chem.MolFromSmiles(other_chain_part)
            if core_mol is None or other_mol is None:
                continue
            try:
                merged = Chem.molzip(core_mol, other_mol)
            except Exception:
                continue

            merged_smiles = Chem.MolToSmiles(merged)
            if merged_smiles.count('[*:') != 1:
                continue
            if '[*:1]' not in merged_smiles:
                merged_smiles = merged_smiles.replace('[*:2]', '[*:1]')

            return {"core": merged_smiles, "target_removed": part}

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
    """규칙의 edit_method에 따라 결합절단형(기존) 또는 원자직접편집형(신규)으로 분기."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    if info.get("edit_method") == "atom_edit":
        from src.tools.atom_editor import apply_atom_edit_from_rule
        return apply_atom_edit_from_rule(smiles, rule_name, candidate_idx)

    located = find_core_and_target(smiles, rule_name)
    if located is None:
        return None
    return reassemble_molecule(located['core'], rule_name, candidate_idx)


def canonicalize(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def iterative_fix_loop(smiles: str, max_iterations: int = 10, candidate_idx: int = 0,
                        llm_client=None, llm_model=None, llm_client_type="gemini"):
    """진단->치환->재평가를 반복.
    llm_client가 주어지면: 어떤 문제부터 고칠지 + 어떤 후보를 쓸지 둘 다 LLM이 판단.
    llm_client_type: "gemini" 또는 "openai_compatible".
    llm_client가 없으면: 리스트 순서(known_problems[0]) + candidate_idx 고정값 사용."""
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
            problem_decision = ask_llm_which_problem_to_fix(llm_client, llm_model, current, problems, client_type=llm_client_type)
            target_rule = problem_decision['rule_name']
            problem_reason = problem_decision.get('reason', '')

            candidate_decision = ask_llm_which_candidate_to_use(llm_client, llm_model, current, target_rule, client_type=llm_client_type)
            chosen_candidate_idx = candidate_decision['candidate_idx']
            candidate_reason = candidate_decision.get('reason', '')
        else:
            target_rule = known_problems[0]['rule_name']
            problem_reason = "규칙 기반(리스트 순서대로)"
            chosen_candidate_idx = candidate_idx
            candidate_reason = "규칙 기반(고정 인덱스)"

        fixed = propose_fix(current, target_rule, chosen_candidate_idx)

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
            "problem_reason": problem_reason,
            "candidate_used": fixed['candidate_used'],
            "candidate_reason": candidate_reason,
        })

    return {"status": "max_iterations_reached", "final_smiles": current, "history": history, "skipped_rules": skipped_rules}
