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
