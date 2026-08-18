"""회귀 방지 테스트 스위트. 반복적으로 겪은 문제들
(함수 중복 정의, 수정이 실제로 반영 안 됨, 핵심 파이프라인 깨짐)을
매 세션 시작 시 한 번에 잡아내기 위한 것. 네트워크/API 호출 없는
테스트만 포함(빠르게, 매번 돌릴 수 있게)."""

import ast
import glob
import os
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _all_source_files():
    patterns = ["src/tools/*.py", "src/models/*.py", "models/*.py"]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(REPO_ROOT, pat)))
    return files


def test_no_duplicate_top_level_function_definitions():
    """%%writefile -a로 같은 함수를 두 번 추가해서, 나중(옛날) 버전이
    최종 반영되는 사고(batch_iterative_fix_loop 사례)를 방지."""
    problems = []
    for path in _all_source_files():
        with open(path) as f:
            source = f.read()
        tree = ast.parse(source)
        names = [node.name for node in tree.body if isinstance(node, ast.FunctionDef)]
        seen = set()
        for name in names:
            if name in seen:
                problems.append(f"{path}: '{name}' 중복 정의")
            seen.add(name)
    assert not problems, "\n".join(problems)


def test_all_source_files_parse():
    """모든 소스 파일이 문법적으로 유효한지."""
    problems = []
    for path in _all_source_files():
        with open(path) as f:
            try:
                ast.parse(f.read())
            except SyntaxError as e:
                problems.append(f"{path}: {e}")
    assert not problems, "\n".join(problems)


def test_precedent_library_structure():
    from src.tools.precedent_library import PRECEDENT_LIBRARY
    assert len(PRECEDENT_LIBRARY) >= 24, f"선례 수가 예상보다 적음: {len(PRECEDENT_LIBRARY)}"
    required_keys = {"rule", "type", "description"}
    for i, p in enumerate(PRECEDENT_LIBRARY):
        missing = required_keys - p.keys()
        assert not missing, f"{i}번째 항목에 키 누락: {missing}"


def test_docking_targets_have_caveat_field():
    """caveat 필드 누락 회귀 방지 (오늘 겪은 EGFR 편향 사고)."""
    from src.tools.docking import DOCKING_TARGETS
    for rule, info in DOCKING_TARGETS.items():
        assert "caveat" in info, f"{rule}에 caveat 필드 없음"


def test_agent_required_functions_exist():
    """agent.py에 있어야 할 핵심 함수/설정 함수들이 다 있는지."""
    import src.tools.agent as agent
    required = [
        "ask_llm_which_problem_to_fix", "ask_llm_which_candidate_to_use",
        "ask_llm_debate_fix", "should_debate",
        "set_debate_budget", "set_sascorer_module", "set_tox_predictor",
        "_try_get_docking_evidence", "_try_compute_score", "_try_get_activity_risk",
    ]
    missing = [name for name in required if not hasattr(agent, name)]
    assert not missing, f"agent.py에 없는 함수: {missing}"


def test_batch_iterative_fix_loop_signature_has_debate_params():
    import inspect
    from src.tools.molecule_editor import batch_iterative_fix_loop
    sig = inspect.signature(batch_iterative_fix_loop)
    assert "use_debate" in sig.parameters
    assert "debate_max_rounds" in sig.parameters


def test_molecule_editor_debate_rounds_recorded_in_history():
    """감사추적용 debate_rounds 키가 iterative_fix_loop 코드에 존재하는지
    (실제 토의 왕복 기록 여부는 별도 통합테스트에서 확인)."""
    with open(os.path.join(REPO_ROOT, "src/tools/molecule_editor.py")) as f:
        content = f.read()
    assert "debate_rounds" in content


@pytest.mark.slow
def test_smoke_iterative_fix_loop_resolves_simple_chain():
    """규칙 기반(LLM 없음) 스모크 테스트: 가장 기본적인 회귀 방지."""
    from src.tools.molecule_editor import iterative_fix_loop, clear_failure_memory
    clear_failure_memory()
    r = iterative_fix_loop("CCCCCCCCCCCCCCCC", max_iterations=10, candidate_idx=0)
    assert r["status"] == "success"


@pytest.mark.slow
def test_smoke_conflict_resolution_picks_lower_remaining_count():
    """다중문제 충돌조정이 실제로 남는 문제 수 적은 쪽을 먼저 고르는지."""
    from src.tools.molecule_editor import iterative_fix_loop, clear_failure_memory
    clear_failure_memory()
    smi = "NNC(=O)CP(=O)(c1ccccc1)c1ccccc1"  # hydrazine + phosphor
    r = iterative_fix_loop(smi, max_iterations=10, candidate_idx=0)
    step1 = r["history"][1]
    assert step1.get("fixed_rule") == "hydrazine"
    assert "충돌 조정" in step1.get("problem_reason", "")
