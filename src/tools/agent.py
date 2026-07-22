import json
from src.tools.replacement_library import get_replacement_candidates


def ask_llm_which_problem_to_fix(client, model_name, smiles, problems):
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

    response = client.models.generate_content(model=model_name, contents=prompt)
    text = response.text.strip()

    if text.startswith('```'):
        text = text.split('```')[1]
        if text.startswith('json'):
            text = text[4:]

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return {"rule_name": known[0]['rule_name'], "reason": "JSON 파싱 실패, 기본값(첫 번째 후보) 사용"}


def ask_llm_which_candidate_to_use(client, model_name, smiles, rule_name):
    """한 문제(rule_name)에 대한 여러 치환 후보 중 어떤 걸 쓸지 LLM에게 판단 요청."""
    info = get_replacement_candidates(rule_name)
    if info is None:
        return None

    candidates = info['candidates']
    if len(candidates) == 1:
        return {"candidate_idx": 0, "reason": "유일한 후보"}

    # LLM에게는 name과 rationale만 보여줌 (smiles는 우리 코드가 나중에 매칭할 참조용)
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

    response = client.models.generate_content(model=model_name, contents=prompt)
    text = response.text.strip()

    if text.startswith('```'):
        text = text.split('```')[1]
        if text.startswith('json'):
            text = text[4:]

    try:
        result = json.loads(text)
        # idx가 유효 범위 안인지 검증
        if not isinstance(result.get('candidate_idx'), int) or not (0 <= result['candidate_idx'] < len(candidates)):
            return {"candidate_idx": 0, "reason": "LLM 응답 idx 범위 오류, 기본값 사용"}
        return result
    except (json.JSONDecodeError, KeyError):
        return {"candidate_idx": 0, "reason": "JSON 파싱 실패, 기본값(첫 번째 후보) 사용"}
