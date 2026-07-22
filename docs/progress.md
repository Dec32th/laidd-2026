# LAIDD 2026 프로젝트 진행 상황
_최종 업데이트: 2026-07-22 (08 노트북 완료 시점)_

## 진행도 요약

| 구성요소 | 상태 | 완료율 |
|---|---|---|
| 환경/인프라 (GitHub-Colab 연동) | 완료 | 100% |
| 도구 계층 (예측/탐지/치환후보/재조립) | 완료 | 100% |
| 반복 루프 (진단→치환→재평가→종료) | 완료 | 100% |
| 에이전트 판단 계층 (문제 우선순위 LLM화) | 완료 | 100% |
| 에이전트 판단 계층 (치환 후보 선택 LLM화) | 미착수 | 0% |
| 치환 라이브러리 확장 (문헌/MMPA 자동) | 미착수 | 0% |
| 평가셋 구성 + 성공률 통계 | 미착수 | 0% |
| 제안서(hwpx) 작성 | 미착수 | 0% |

**전체 진행도: 약 55%** (8개 항목 중 4개 완료, 4개 미착수 기준 단순 평균)

## 🎯 최종 목표
JUMP AI 2026 예선 (8/7 마감) — "도구 활용 기반 분자 최적화 루프" 분야
설계: 도구 계층(결정적) + 에이전트 계층(LLM 판단)의 하이브리드 구조로,
독성 있는 분자를 진단→치환 후보 판단→재조립→재평가하는 루프

## ✅ 완료된 것

### 환경/인프라
- GitHub repo(Dec32th/laidd-2026) + Colab 연동 (private repo, PAT 토큰 인증)
- 폴더 구조: src/tools/, notebooks/(번호순), models/, docs/

### 도구 계층
| 도구 | 파일 | 역할 | 결과 |
|---|---|---|---|
| 데이터 전처리 | data_prep.py | Tox21 로드+필터링(8개 제외)+SMILES 반환 | 완료 |
| 독성 예측 baseline | (모델) | ECFP+RandomForest, 12개 assay | ==평균 AUROC 0.821== |
| ChemBERTa 비교 | 03 노트북 | 사전학습 임베딩 비교실험 | AUROC 0.787, ==baseline 채택을 실험적으로 정당화한 근거==🔑 |
| 문제구조 탐지 | toxicophore_detector.py | FilterCatalog(PAINS+BRENK) → rule_name+원자인덱스 | 완료 |
| 치환 후보 라이브러리 | replacement_library.py | 7개 규칙(nitro, aldehyde, michael_acceptor, thiourea, acyl_halide, alkyl_halide, aniline) → 각 후보+rationale | 완료 |
| 분자 재조립 | molecule_editor.py | find_core_and_target + reassemble_molecule + propose_fix | 완료, ==실제 화학적으로 유효한 신규 분자 생성 검증됨==(예: 클로로아세트산→글리콜산) |

### 반복 루프 + 에이전트 판단 계층
- `iterative_fix_loop()`: 진단→치환→재평가 반복, 종료조건 4종
  (success / stuck / cycle_detected / no_known_fix / max_iterations_reached)
- ==canonicalize()로 SMILES 정규화 후 순환 감지 — "같은 분자, 다른 표기"로 인한 오탐 방지==
- ==모르는 규칙은 멈추지 않고 건너뛰며 skipped_rules에 투명하게 기록==🔑
  (연구윤리/투명성 평가 항목과 직결)
- `agent.py`: 여러 toxicophore 중 어떤 것부터 고칠지 **LLM(Gemini 3.5 Flash)이 판단**
  - 후보가 1개뿐이면 LLM 호출 생략 (효율성)
  - JSON 강제 출력 + 파싱 실패시 규칙기반 fallback으로 안전장치 확보
- ==**규칙기반 vs LLM기반 비교실험 완료**==🔑: 같은 분자(aniline+nitro_group 동반 보유)에서
  - 규칙기반: aniline부터 고쳐 1스텝 만에 success
  - LLM기반: "니트로기가 아닐린보다 유전독성 위험이 크다"는 화학적 근거로 nitro_group을 먼저 선택, 2스텝 거쳐 success
  - → ==단순 순서가 아닌 화학적 판단이 실제로 다른 경로를 만든다는 것을 실증==

## 🔲 남은 것 (우선순위 순)
1. 치환 후보(candidate) 선택도 LLM이 판단하도록 확장 (현재는 candidate_idx=0 고정)
2. 치환 라이브러리 확장 — 두 갈래로 고려 중:
   - Medicinal chemistry 리뷰논문/사례 기반 수동 추가 (화학적 신뢰도 높음, 항목당 시간 소요)
   - ==rdMMPA를 활용해 Tox21 데이터에서 자동으로 치환쌍(matched molecular pair) 추출==🔑
     (데이터 기반 자동 확장 — "확장 가능한 구조"로 제안서에 서술 가능한 포인트)
3. Held-out set 전체에 루프 적용 → 성공률/개선율 통계
4. 실제 사례(개발 중단 약물 등) 케이스스터디 1~2개 준비
5. 제안서(hwpx) 작성

## ⚠️ 시스템 한계 (docs/limitations.md 요약)
- ==정적 구조 분석만 수행 — 실제 대사(metabolism) 시뮬레이션은 하지 않음==
- 약물-약물 상호작용은 스코프 밖 (DrugBank 등 별도 데이터셋 필요)
- 조합 효과는 사후 재검사로만 대응 (사전 예측 안 함)
- ==본 도구는 hit-to-lead 초기 스크리닝 보조 도구로 위치 설정 — 최종 검증엔 ADMET/실험 필요==
  (→ 한계를 명확히 인지하고 서술한 것 자체가 완성도 평가에 플러스 요인)

## ⚠️ 자주 헷갈렸던 것 (작업 습관 메모)
- 새 Colab 세션마다: pip install → git clone(토큰) → %cd 경로 확인 → git config 재설정 필요
- 파일 저장 전엔 항상 !pwd로 위치 확인 습관화
- 모듈 수정 후엔 importlib.reload() 필요 (안 하면 옛날 버전 계속 씀)
- SMARTS 패턴 설계 시 "패턴 포함 여부"와 "정확한 크기 매칭"은 다름 (alkyl_halide 버그 사례)