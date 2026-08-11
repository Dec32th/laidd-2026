from rdkit import Chem
from rdkit.Chem import rdMMPA
from src.tools.replacement_library import get_replacement_candidates
import hashlib

def _library_version_hash():
    from src.tools.replacement_library import get_replacement_candidates
    lib = get_replacement_candidates.__globals__['REPLACEMENT_LIBRARY']
    content_str = str(sorted(lib.items()))
    return hashlib.md5(content_str.encode()).hexdigest()[:8]

_FAILURE_MEMORY = {}


def clear_failure_memory():
    global _FAILURE_MEMORY
    _FAILURE_MEMORY = {}


def _check_and_match(part_smiles, problem_pattern, pattern_size):
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


def _candidate_order_for_rule(rule_name: str, preferred_idx: int):
    info = get_replacement_candidates(rule_name)
    if info is None:
        return [preferred_idx]
    n = len(info['candidates'])
    order = [preferred_idx] if 0 <= preferred_idx < n else []
    order += [i for i in range(n) if i != preferred_idx]
    return order


def iterative_fix_loop(smiles: str, max_iterations: int = 10, candidate_idx: int = 0,
                        llm_client=None, llm_model=None, llm_client_type="gemini",
                        use_failure_memory: bool = True):
    """진단->치환->재평가를 반복.
    핵심: candidate가 '화학적으로 유효(is_valid)'해도 대상 규칙이 실제로
    해소됐는지 재진단(detect_toxicophores)까지 확인한다. 그렇지 않으면
    항상 valid하지만 문제를 안 고치는 candidate(예: 단순 삽입형)가
    무한 반복 채택되어 진짜 해법(예: 분기형)으로 넘어가지 못하는 문제가
    있었음. 완전 해소가 안 되면 마지막으로 시도한(=대개 더 나은)
    valid 결과를 fallback으로 채택해 다음 iteration에서 계속 개선."""
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
                preferred_candidate_idx = candidate_decision['candidate_idx']
                this_candidate_reason = candidate_decision.get('reason', '')

                if preferred_candidate_idx == -1:
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
                preferred_candidate_idx = candidate_idx
                this_candidate_reason = "규칙 기반(고정 인덱스 우선, 실패/미해소 시 같은 규칙 내 다른 candidate로 재시도)"

            rule_fixed = None
            fallback_attempt = None
            fallback_used_idx = None
            fallback_reason = None

            for try_idx in _candidate_order_for_rule(candidate_rule, preferred_candidate_idx):
                memory_key = (current, candidate_rule, try_idx, _library_version_hash())
                if use_failure_memory and memory_key in _FAILURE_MEMORY:
                    failed_attempts.append(f"{candidate_rule}[idx={try_idx}](memory-skip)")
                    continue

                attempt = propose_fix(current, candidate_rule, try_idx)
                if attempt is None or not attempt.get('is_valid'):
                    failed_attempts.append(f"{candidate_rule}[idx={try_idx}]")
                    if use_failure_memory:
                        _FAILURE_MEMORY[memory_key] = True
                    continue

                # valid해도 실제로 이 규칙이 재진단에서 사라졌는지 확인
                recheck = detect_toxicophores(attempt['new_smiles'])
                still_flagged = any(p['rule_name'] == candidate_rule for p in recheck)

                if not still_flagged:
                    rule_fixed = attempt
                    candidate_reason = f"{this_candidate_reason} (candidate_idx={try_idx}, 완전 해소)"
                    break
                else:
                    failed_attempts.append(f"{candidate_rule}[idx={try_idx}](valid이나 미해소)")
                    fallback_attempt = attempt
                    fallback_used_idx = try_idx
                    fallback_reason = f"{this_candidate_reason} (candidate_idx={try_idx}, 부분 개선/다음 iteration에서 계속)"

            if rule_fixed is None and fallback_attempt is not None:
                rule_fixed = fallback_attempt
                candidate_reason = fallback_reason

            if rule_fixed is not None:
                fixed = rule_fixed
                target_rule = candidate_rule
                break

        if fixed is None:
            reason_detail = (f"이 단계에서 known 규칙들의 모든 candidate를 순서대로 시도했으나 "
                              f"({failed_attempts}) 모두 실행에 실패했습니다(memory-skip 표시는 이전에 "
                              f"실패했던 것으로 확인되어 재시도 없이 건너뛴 항목). 흔한 원인: 유기금속/무기염 "
                              f"등 특수 화학종, 고리 구조와의 예상치 못한 충돌, 또는 원자가 계산 오류입니다.")
            return {"status": "stuck", "reason": f"시도한 규칙/candidate {failed_attempts} 모두 치환 실패",
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


def batch_iterative_fix_loop(smiles_list, max_iterations=10, candidate_idx=0,
                               llm_client=None, llm_model=None, llm_client_type="gemini",
                               max_workers=5, progress=True):
    """여러 분자에 iterative_fix_loop를 스레드 병렬로 적용.
    LLM API 호출이 병목인 경우(네트워크 대기 시간) 유효한 개선이며,
    화학 계산 로직(iterative_fix_loop 자체)은 전혀 수정하지 않는다.
    반환: [(smiles, result_dict), ...] (완료 순서, 입력 순서와 다를 수 있음)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _process_one(smi):
        r = iterative_fix_loop(
            smi, max_iterations=max_iterations, candidate_idx=candidate_idx,
            llm_client=llm_client, llm_model=llm_model, llm_client_type=llm_client_type,
        )
        return smi, r

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_one, smi): smi for smi in smiles_list}
        for i, future in enumerate(as_completed(futures)):
            smi, r = future.result()
            results.append((smi, r))
            if progress:
                print(f"[{i+1}/{len(smiles_list)}] {smi[:30]} -> {r['status']}")
    return results


def batch_iterative_fix_loop(smiles_list, max_iterations=10, candidate_idx=0,
                               llm_client=None, llm_model=None, llm_client_type="gemini",
                               max_workers=5, progress=True):
    """여러 분자에 iterative_fix_loop를 스레드 병렬로 적용.
    LLM API 호출이 병목인 경우(네트워크 대기 시간) 유효한 개선이며,
    화학 계산 로직(iterative_fix_loop 자체)은 전혀 수정하지 않는다.
    반환: [(smiles, result_dict), ...] (완료 순서, 입력 순서와 다를 수 있음)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _process_one(smi):
        r = iterative_fix_loop(
            smi, max_iterations=max_iterations, candidate_idx=candidate_idx,
            llm_client=llm_client, llm_model=llm_model, llm_client_type=llm_client_type,
        )
        return smi, r

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_one, smi): smi for smi in smiles_list}
        for i, future in enumerate(as_completed(futures)):
            smi, r = future.result()
            results.append((smi, r))
            if progress:
                print(f"[{i+1}/{len(smiles_list)}] {smi[:30]} -> {r['status']}")
    return results
