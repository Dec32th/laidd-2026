import json
import time
from src.tools.replacement_library import get_replacement_candidates

_llm_error_log = []
_llm_consecutive_failures = 0
_LLM_FAILURE_LIMIT = 3
_debate_call_budget = {"remaining": 100}

def set_debate_budget(n):
    """토의(debate)에 쓸 수 있는 총 LLM 호출 수 상한을 재설정."""
    _debate_call_budget["remaining"] = n

_injected_sascorer = {"module": None}
_injected_tox_predictor = {"fn": None}

def set_sascorer_module(module):
    """세션마다 다운로드/import한 sascorer 모듈을 등록."""
    _injected_sascorer["module"] = module

def set_tox_predictor(fn):
    """(original_smiles, fixed_smiles, rule_name) -> tox_delta(float) 또는 None
    을 반환하는 콜백을 등록. 노트북마다 다르게 학습한 baseline 모델을
    감싸서 넘기면 됨."""
    _injected_tox_predictor["fn"] = fn

def _call_llm(client, model_name, prompt, client_type="gemini"):
    """client_type에 따라 Gemini SDK 또는 OpenAI 호환 SDK로 호출하고,
    응답 텍스트만 통일된 형태로 반환."""
    if client_type == "gemini":
        response = client.models.generate_content(model=model_name, contents=prompt)
        return response.text
    elif client_type == "openai_compatible":
        for attempt in range(2):  # rate limit 시 1회만 재시도
            try:
                response = client.chat.completions.create(
                    model=model_name,
                    messages=[{"role": "user", "content": prompt}],
                    max_tokens=500,
                    timeout=30,
                    extra_body={"enable_thinking": False},
                )
                return response.choices[0].message.content
            except Exception as e:
                _llm_error_log.append(repr(e))
                if 'RateLimitError' in type(e).__name__ and attempt == 0:
                    time.sleep(3)
                    continue
                return f"ERROR: LLM 호출 실패/타임아웃 - {e}"
    else:
        raise ValueError(f"알 수 없는 client_type: {client_type}")


def _parse_json_response(text, fallback):
    """LLM 응답에서 JSON을 추출/파싱. 마크다운 코드블록(```json)으로
    감싸져 온 경우를 벗겨내고, 파싱 실패 시 fallback을 반환."""
    text = text.strip()
    if text.startswith('```'):
        text = text.split('```')[1]
        if text.startswith('json'):
            text = text[4:]
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return fallback

def _try_get_docking_evidence(rule_name, smiles_before, smiles_after):
    """도킹 표적이 등록된 규칙이면 자동으로 도킹 실행, 아니면 None.
    도킹 실패/미등록/예외는 전부 조용히 None으로 처리(critic 프롬프트에서
    도킹 근거 없이 진행하는 것으로 자연스럽게 폴백)."""
    try:
        from src.tools.docking import auto_dock_precedent, DOCKING_TARGETS
        if rule_name not in DOCKING_TARGETS:
            return None
        result = auto_dock_precedent(rule_name, smiles_before, smiles_after)
        if result.get('error') or result.get('delta') is None:
            return None
        return result
    except Exception:
        return None

def _try_get_warhead_reference(rule_name):
    """워헤드 반응성 참고표가 있으면 조회, 없으면 조용히 None."""
    try:
        from src.tools.docking import get_warhead_reference
        return get_warhead_reference(rule_name)
    except Exception:
        return None

def _try_compute_score(rule_name, smiles_before, smiles_after, docking_evidence=None):
    """등록된 sascorer/tox_predictor/도킹 결과를 모아 종합 점수 계산.
    일부만 등록돼 있어도 compute_multi_objective_score가 나머지로
    자동 정규화하므로 실패하지 않음. 계산 자체가 실패하면 None."""
    try:
        from src.tools.scoring import compute_multi_objective_score

        tox_delta = None
        if _injected_tox_predictor["fn"] is not None:
            try:
                tox_delta = _injected_tox_predictor["fn"](smiles_before, smiles_after, rule_name)
            except Exception:
                tox_delta = None

        precedent_docking_delta = docking_evidence["delta"] if docking_evidence else None

        return compute_multi_objective_score(
            smiles_before, smiles_after, rule_name,
            tox_delta=tox_delta,
            sascorer_module=_injected_sascorer["module"],
            precedent_docking_delta=precedent_docking_delta,
        )
    except Exception:
        return None

def _try_get_activity_risk(rule_name, smiles_before, smiles_after):
    """활성 보존 위험도 평가. sascorer가 등록 안 돼 있으면 조용히 None
    (activity_metrics가 sascorer_module을 필수로 요구하므로)."""
    try:
        if _injected_sascorer["module"] is None:
            return None
        from src.tools.activity_metrics import (compute_activity_preservation_metrics,
                                                   classify_activity_risk_v3)
        metrics = compute_activity_preservation_metrics(
            smiles_before, smiles_after, _injected_sascorer["module"]
        )
        if metrics is None:
            return None
        return classify_activity_risk_v3(metrics)
    except Exception:
        return None

def ask_llm_which_problem_to_fix(client, model_name, smiles, problems, client_type="gemini"):
    """여러 toxicophore 중 어떤 것부터 고칠지 LLM에게 판단을 요청."""
    known = [p for p in problems if get_replacement_candidates(p['rule_name']) is not None]

    if not known:
        return None
    if len(known) == 1:
        return {"rule_name": known[0]['rule_name'], "reason": "유일한 치환 가능 후보"}

    prompt = f"""당신은 신약개발 화학자입니다. 다음 분자에서 여러 구조적 문제(toxicophore)가 발견되었습니다.

분자 SMILES: {smiles}

발견된 문제 중, 우리가 실제로 치환 가능한 것들:
{json.dumps(known, ensure_ascii=False, indent=2)}

이 중 어떤 문제를 먼저 해결하는 것이 화학적으로 더 타당한지 판단하고,
반드시 아래 JSON 형식으로만 답하세요. 다른 설명 없이 JSON만 출력하세요.

{{"rule_name": "선택한 문제의 rule_name", "reason": "선택 이유 한 문장"}}
"""

    text = _call_llm(client, model_name, prompt, client_type)
    fallback = {"rule_name": known[0]['rule_name'], "reason": "JSON 파싱 실패, 기본값(첫 번째 후보) 사용"}
    return _parse_json_response(text, fallback)


def ask_llm_which_candidate_to_use(client, model_name, smiles, rule_name, client_type="gemini"):
    """한 문제(rule_name)에 대한 여러 치환 후보 중 어떤 걸 쓸지 LLM에게 판단 요청.

    candidate의 rationale 중 하나라도 '[참고]'로 시작하는 문구가 있으면,
    이는 실제 승인약물 사례에서 이 골격이 안전하게 쓰인 경우가 있다는 뜻이므로,
    candidate가 1개뿐이더라도(원래는 LLM 호출을 건너뛰던 경우) 반드시 LLM에게
    판단을 맡긴다. 이 경우 LLM은 candidate_idx로 -1을 반환하여 "치환을
    보류하고 사람(연구자) 검토가 필요하다"고 명시적으로 표시할 수 있다.
    """
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    candidates = info['candidates']
    has_caution = any('[참고]' in c.get('rationale', '') for c in candidates)

    if len(candidates) == 1 and not has_caution:
        return {"candidate_idx": 0, "reason": "유일한 후보"}

    candidate_info = [
        {"idx": i, "name": c['name'], "rationale": c['rationale']}
        for i, c in enumerate(candidates)
    ]

    prompt = f"""당신은 신약개발 화학자입니다. 다음 분자에서 '{rule_name}' 문제를
해결하기 위한 치환 후보가 있습니다.

분자 SMILES: {smiles}

치환 후보들:
{json.dumps(candidate_info, ensure_ascii=False, indent=2)}

각 후보의 rationale에 "[참고]"로 시작하는 문구가 있다면, 이는 "이 골격이
실제 승인 약물에서 반응성이 아닌 안정적 형태로 널리 쓰인 사례가 있으니,
경고를 절대적 기준이 아닌 참고 신호로 해석하라"는 뜻입니다. 이 경우 먼저
"이 분자가 그 참고사항이 가리키는 안전한 사용 사례와 실제로 유사한지"를
판단하세요.
- 유사하다고 판단되면서, 후보가 여러 개라면 변화 폭이 더 작은 후보를 선택하세요.
- 유사하다고 판단되고, 치환 자체가 불필요하다고 볼 만큼 뚜렷하다면,
  candidate_idx를 -1로 답해 "치환 보류, 사람 검토 필요"를 표시하세요.
- 참고사항이 없거나 이 분자가 그 사례와 유사하지 않다면, 평소대로 가장
  적절한 후보를 선택하세요.

반드시 아래 JSON 형식으로만 답하세요. 다른 설명 없이 JSON만 출력하세요.

{{"candidate_idx": 선택한 후보의 idx(정수, 또는 보류 시 -1), "reason": "판단 이유 한 문장"}}
"""

    text = _call_llm(client, model_name, prompt, client_type)
    fallback = {"candidate_idx": 0, "reason": "JSON 파싱 실패, 기본값(첫 번째 후보) 사용"}
    result = _parse_json_response(text, fallback)

    idx = result.get('candidate_idx')
    if not isinstance(idx, int) or not (-1 <= idx < len(candidates)):
        return {"candidate_idx": 0, "reason": "LLM 응답 idx 범위 오류, 기본값 사용"}
    return result


def ask_llm_debate_fix(client, model_name, smiles_before, smiles_after, rule_name,
                        candidate_name, candidate_rationale, client_type="gemini",
                        max_rounds=2):
    """제안자(원래 candidate를 고른 논리)와 검토자(critic)가 여러 라운드
    대화하며 합의에 도달하려 시도. 매 라운드 critic이 판단하고, 반려하면
    proposer가 반박, critic이 재판단. max_rounds 안에 합의(양쪽 다 승인,
    또는 critic이 최종 반려로 확정) 안 되면 "escalate"로 사람 검토行.

    반환: {"final_verdict": "approved"|"rejected"|"escalate",
           "rounds": [{"role": "critic"|"proposer", "text": str}, ...],
           "consensus_reached": bool}
    """
    if _debate_call_budget["remaining"] <= 0:
        return {"final_verdict": "approved", "rounds": [], "consensus_reached": True,
                "budget_exhausted": True}
    _debate_call_budget["remaining"] -= 1
    rounds_log = []
    proposer_argument = candidate_rationale

    for round_num in range(1, max_rounds + 1):
        docking_evidence = _try_get_docking_evidence(rule_name, smiles_before, smiles_after)
        score_result = _try_compute_score(rule_name, smiles_before, smiles_after, docking_evidence)
        activity_risk = _try_get_activity_risk(rule_name, smiles_before, smiles_after)
        warhead_ref = _try_get_warhead_reference(rule_name)
        critic_prompt = f"""당신은 신약개발 화학 검토자(critic)입니다. 동료 화학자가 아래
치환을 제안했습니다.

원본 분자: {smiles_before}
치환 후 분자: {smiles_after}
해결하려던 문제: {rule_name}
제안된 치환: {candidate_name}
제안자의 근거: {proposer_argument}
{f"실측 도킹 결합력 변화: {docking_evidence['target']} 표적, {docking_evidence['score_original']:.2f} → {docking_evidence['score_fixed']:.2f} kcal/mol (delta {docking_evidence['delta']:+.2f}). 이 정량 데이터를 판단에 반영하세요.{' [주의: ' + docking_evidence['caveat'] + ']' if docking_evidence and docking_evidence.get('caveat') else ''}" if docking_evidence else ""}
{f"종합 점수: {score_result['composite_score']:.2f} (세부: {score_result['component_scores']}). 이것도 판단에 참고하세요." if score_result and score_result.get('composite_score') is not None else ""}
{f"활성 보존 위험도 평가: {activity_risk['verdict']} (세부: {'; '.join(activity_risk['details'])})" if activity_risk else ""}
{f"워헤드 반응성 참고: {warhead_ref['reactivity_note']} 실무 지침: {warhead_ref['practical_guidance']}" if warhead_ref else ""}

이 치환에 동의하는지 비판적으로 검토하세요. 동의하지 않는다면 구체적으로
어떤 점이 문제인지 명시하세요(새로운 독성 구조 생성 가능성, 근거의
논리적 결함, precedent 오독 등).

반드시 아래 JSON 형식으로만 답하세요.
{{"verdict": "approved" 또는 "rejected", "reason": "판단 이유, 반려 시 구체적 반론 포함"}}
"""
        critic_text = _call_llm(client, model_name, critic_prompt, client_type)
        critic_result = _parse_json_response(
            critic_text, {"verdict": "approved", "reason": "JSON 파싱 실패, 기본 승인"}
        )
        rounds_log.append({"role": "critic", "round": round_num, "text": critic_result})

        if critic_result.get("verdict") == "approved":
            return {"final_verdict": "approved", "rounds": rounds_log, "consensus_reached": True}

        if round_num == max_rounds:
            break

        proposer_prompt = f"""당신은 방금 아래 치환을 제안한 화학자입니다.

원본 분자: {smiles_before}
치환 후 분자: {smiles_after}
당신의 원래 근거: {proposer_argument}

동료 검토자(critic)가 다음과 같이 반려했습니다: "{critic_result.get('reason', '')}"

이 반론에 대해 답하세요. 반론이 타당하면 인정하고 제안을 철회하세요.
반론이 부당하다면 왜 원래 치환이 여전히 타당한지 반박하세요.

반드시 아래 JSON 형식으로만 답하세요.
{{"stance": "withdraw" 또는 "defend", "argument": "반박 또는 철회 이유"}}
"""
        proposer_text = _call_llm(client, model_name, proposer_prompt, client_type)
        proposer_result = _parse_json_response(
            proposer_text, {"stance": "withdraw", "argument": "JSON 파싱 실패, 기본 철회"}
        )
        rounds_log.append({"role": "proposer", "round": round_num, "text": proposer_result})

        if proposer_result.get("stance") == "withdraw":
            return {"final_verdict": "rejected", "rounds": rounds_log, "consensus_reached": True}

        proposer_argument = proposer_result.get("argument", proposer_argument)

    return {"final_verdict": "escalate", "rounds": rounds_log, "consensus_reached": False}
def should_debate(candidate_rationale):
    """이 candidate가 토의(debate)를 거칠 필요가 있는지 판단.
    rationale에 '[참고]'가 있으면 실제 승인약물 사례와 겹칠 수 있다는
    뜻이므로, 단순 채택 대신 토의로 한 번 더 검토해야 함."""
    return '[참고]' in (candidate_rationale or '')

def ask_llm_debate_fix_consistent(client, model_name, smiles_before, smiles_after, rule_name,
                                    candidate_name, candidate_rationale, client_type="gemini",
                                    max_rounds=2, n_repeats=3):
    """ask_llm_debate_fix를 n_repeats번 반복해서 다수결로 최종 판정.
    LLM 토의의 확률적 특성(같은 입력에도 매번 다른 판정)을 완화하기 위함.
    반환: {"final_verdict": ..., "vote_counts": {...}, "all_results": [...]}
    """
    from collections import Counter
    results = []
    for _ in range(n_repeats):
        r = ask_llm_debate_fix(
            client, model_name, smiles_before, smiles_after, rule_name,
            candidate_name, candidate_rationale, client_type=client_type, max_rounds=max_rounds,
        )
        results.append(r)

    votes = Counter(r['final_verdict'] for r in results)
    majority_verdict, majority_count = votes.most_common(1)[0]

    return {
        "final_verdict": majority_verdict,
        "vote_counts": dict(votes),
        "consensus_strength": f"{majority_count}/{n_repeats}",
        "all_results": results,
    }
