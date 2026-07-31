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
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    problem_pattern = Chem.MolFromSmarts(info['problem_smarts'])
    pattern_size = problem_pattern.GetNumAtoms()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

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
    """진단->치환->재평가를 반복. known 규칙 중 우선순위가 가장 높은 것이
    propose_fix에서 실패하면, stuck 처리 전에 같은 분자의 다른 known
    규칙들을 순서대로 시도한다."""
    from src.tools.toxicophore_detector import detect_toxicophores
    from src.tools.agent import ask_llm_which_problem_to_fix, ask_llm_which_candidate_to_use

    current = canonicalize(smiles)
    seen = {current}
    history = [{"step": 0, "smiles": current}]
    skipped_rules = []
    skipped_details = []
    flagged_for_review = set()

    for step in range(1, max_iterations + 1):
        problems = detect_toxicophores(current)
        history[-1]["problems"] = problems

        if not problems:
            return {"status": "success", "final_smiles": current, "history": history,
                    "skipped_rules": skipped_rules, "skipped_details": skipped_details}

        known_problems = [p for p in problems if get_replacement_candidates(p['rule_name']) is not None
                          and p['rule_name'] not in flagged_for_review]
        unknown_problems = [p for p in problems if get_replacement_candidates(p['rule_name']) is None]

        for p in unknown_problems:
            if p['rule_name'] not in skipped_rules:
                skipped_rules.append(p['rule_name'])
                mol_cur = Chem.MolFromSmiles(current)
                matched_atoms = p['atom_indices']
                atom_symbols = [mol_cur.GetAtomWithIdx(i).GetSymbol() for i in matched_atoms] if mol_cur else []
                skipped_details.append({
                    "rule_name": p['rule_name'],
                    "reason": f"라이브러리에 등록되지 않은 규칙입니다. FilterCatalog(PAINS/BRENK)가 "
                              f"'{p['rule_name']}'로 진단했으며, 매치된 원자 인덱스는 {matched_atoms}"
                              f"(원소: {atom_symbols})입니다. 이 구조에 대한 치환 규칙을 "
                              f"replacement_library.py에 추가하면 자동으로 처리 가능합니다.",
                    "atom_indices": matched_atoms,
                })

        if not known_problems:
            return {"status": "no_known_fix", "final_smiles": current, "history": history,
                    "skipped_rules": skipped_rules, "skipped_details": skipped_details}

        if llm_client is not None:
            problem_decision = ask_llm_which_problem_to_fix(llm_client, llm_model, current, problems, client_type=llm_client_type)
            preferred_rule = problem_decision['rule_name']
            problem_reason = problem_decision.get('reason', '')
            ordered_rules = [preferred_rule] + [p['rule_name'] for p in known_problems if p['rule_name'] != preferred_rule]
        else:
            problem_reason = "규칙 기반(리스트 순서대로)"
            ordered_rules = [p['rule_name'] for p in known_problems]

        fixed = None
        target_rule = None
        candidate_reason = None
        failed_attempts = []

        for candidate_rule in ordered_rules:
            if llm_client is not None:
                candidate_decision = ask_llm_which_candidate_to_use(llm_client, llm_model, current, candidate_rule, client_type=llm_client_type)
                chosen_candidate_idx = candidate_decision['candidate_idx']
                this_candidate_reason = candidate_decision.get('reason', '')

                if chosen_candidate_idx == -1:
                    flagged_for_review.add(candidate_rule)
                    if candidate_rule not in skipped_rules:
                        skipped_rules.append(candidate_rule)
                    skipped_details.append({
                        "rule_name": candidate_rule,
                        "reason": f"LLM이 치환을 보류했습니다: {this_candidate_reason} "
                                  f"(이 분자가 [참고] 사항에 해당하는 안전한 실사용 사례와 유사하다고 "
                                  f"판단되어, 자동 치환 대신 연구자의 직접 검토를 권장합니다.)",
                        "atom_indices": next((p['atom_indices'] for p in problems if p['rule_name'] == candidate_rule), []),
                    })
                    continue
            else:
                chosen_candidate_idx = candidate_idx
                this_candidate_reason = "규칙 기반(고정 인덱스)"

            attempt = propose_fix(current, candidate_rule, chosen_candidate_idx)
            if attempt is not None and attempt.get('is_valid'):
                fixed = attempt
                target_rule = candidate_rule
                candidate_reason = this_candidate_reason
                break
            else:
                failed_attempts.append(candidate_rule)

        if fixed is None:
            reason_detail = (f"이 단계에서 known 규칙 {failed_attempts} 전부를 순서대로 시도했으나 "
                              f"모두 실행에 실패했습니다. 흔한 원인: 유기금속/무기염 등 특수 화학종, "
                              f"고리 구조와의 예상치 못한 충돌, 또는 원자가 계산 오류입니다.")
            return {"status": "stuck", "reason": f"시도한 규칙 {failed_attempts} 모두 치환 실패",
                    "reason_detail": reason_detail,
                    "final_smiles": current, "history": history,
                    "skipped_rules": skipped_rules, "skipped_details": skipped_details}

        new_current = canonicalize(fixed['new_smiles'])

        if new_current in seen:
            return {"status": "cycle_detected", "final_smiles": current, "history": history,
                    "skipped_rules": skipped_rules, "skipped_details": skipped_details}

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

    return {"status": "max_iterations_reached", "final_smiles": current, "history": history,
            "skipped_rules": skipped_rules, "skipped_details": skipped_details}
