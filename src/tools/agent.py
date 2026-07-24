import json
from src.tools.replacement_library import get_replacement_candidates


def _call_llm(client, model_name, prompt, client_type="gemini"):
    """client_type에 따라 Gemini SDK 또는 OpenAI 호환 SDK로 호출하고,
    응답 텍스트만 통일된 형태로 반환."""
    if client_type == "gemini":
        response = client.models.generate_content(model=model_name, contents=prompt)
        return response.text
    elif client_type == "openai_compatible":
        response = client.chat.completions.create(
            model=model_name,
            messages=[{"role": "user", "content": prompt}]
        )
        return response.choices[0].message.content
    else:
        raise ValueError(f"알 수 없는 client_type: {client_type}")


def _parse_json_response(text, fallback):
    text = text.strip()
    if text.startswith('```'):
        text = text.split('```')[1]
        if text.startswith('json'):
            text = text[4:]
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return fallback


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
    """한 문제(rule_name)에 대한 여러 치환 후보 중 어떤 걸 쓸지 LLM에게 판단 요청."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    candidates = info['candidates']
    if len(candidates) == 1:
        return {"candidate_idx": 0, "reason": "유일한 후보"}

    candidate_info = [
        {"idx": i, "name": c['name'], "rationale": c['rationale']}
        for i, c in enumerate(candidates)
    ]

    prompt = f"""당신은 신약개발 화학자입니다. 다음 분자에서 '{rule_name}' 문제를
해결하기 위한 여러 치환 후보가 있습니다.

분자 SMILES: {smiles}

치환 후보들:
{json.dumps(candidate_info, ensure_ascii=False, indent=2)}

이 중 이 분자 맥락에서 가장 적절한 후보를 선택하고,
반드시 아래 JSON 형식으로만 답하세요. 다른 설명 없이 JSON만 출력하세요.

{{"candidate_idx": 선택한 후보의 idx(정수), "reason": "선택 이유 한 문장"}}
"""

    text = _call_llm(client, model_name, prompt, client_type)
    fallback = {"candidate_idx": 0, "reason": "JSON 파싱 실패, 기본값(첫 번째 후보) 사용"}
    result = _parse_json_response(text, fallback)

    if not isinstance(result.get('candidate_idx'), int) or not (0 <= result['candidate_idx'] < len(candidates)):
        return {"candidate_idx": 0, "reason": "LLM 응답 idx 범위 오류, 기본값 사용"}
    return result
