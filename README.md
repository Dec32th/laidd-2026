# ToxGuard — Toxicophore-Guided Molecular Substitution System

Tox21 데이터셋의 분자에서 독성 구조(toxicophore)를 자동 진단하고,
화학적으로 타당한 치환을 제안·검증하는 파이프라인. JUMP AI 2026
공모전(4th JUMP AI competition, 팀명 MorForge) 제출작을 기반으로
지속 개발 중.

## 개요

"LLM을 tool-calling 에이전트로 활용해 도메인 규칙 기반 진단과 결합한 시스템"
RDKit 규칙 기반 진단 + LLM 기반 토의(proposer-critic) 검증 + 단백질
도킹/합성용이성 정량 근거를 결합한 다단계 시스템:

1. 진단: BRENK/PAINS 필터로 분자 내 toxicophore 탐지
2. 치환 후보 생성: 규칙별 치환 라이브러리(replacement_library.py)
   + 다중문제 충돌조정(그리디, 남는 문제 수 최소화)
3. 토의 검증: candidate에 [참고] 표시(승인약물 사례와 유사할
   가능성)가 있으면 LLM proposer-critic 다라운드 토의로 재검토
4. 정량 근거: 단백질 도킹(AutoDock Vina), 합성용이성(SA score),
   활성보존 지표(Tanimoto, QED, LogP 등)를 critic 판단에 제공
5. 감사추적: 모든 판단(승인/반려/에스컬레이트)과 근거를 기록

## 최종 수치 (valid set 1173개)

| | 완전해결 |
|---|---|
| 예선 제출 | 52.0% |
| 규칙 기반 (현재) | 74.5% |
| 토의 파이프라인 (현재) | 67~70% |

토의 파이프라인의 수치가 더 낮은 건 성능 저하가 아니라 의도된
트레이드오프 — 화학적으로 의심스러운 자동 치환을 막아 표면적
성공률은 낮아지지만 각 성공의 신뢰도가 올라감.

## 주요 발견 사례

- BRENK SMARTS 실측 기반 재설계: Aliphatic_long_chain 규칙이
  원소 종류와 무관한 순수 위상학적 패턴임을 직접 분석으로 확인,
  전략 전환으로 완전해결
- 도킹 편향 발견 및 수정: 표적 특이성 없는 벤치마크 도킹 결과를
  critic이 이 분자는 이 표적이다로 오귀속하는 시스템적 편향 발견,
  caveat 메커니즘으로 수정
- 정직한 실패와 재도전: 표준 도킹의 공유결합 형성 미반영 문제를
  거리 기반으로 근사하려 두 차례 시도해 실패, 이후 Meeko+AutoDock-GPU
  기반으로 근본 해결 성공(오시메르티닙 -5.31 vs 이타콘산 +7.06
  kcal/mol, docs/meeko_covalent_docking_progress.md 참고)
- 파괴적 편집 가드: 다중문제 충돌조정이 분자를 통째로 파괴하는
  candidate를 최선으로 착각할 수 있는 위험을 데이터로 발견하고 수정

더 자세한 개발 과정은 CHANGELOG.md와 docs/limitations.md를 참고.

## 저장소 구조

- src/tools/ : 핵심 로직 (진단, 치환, 토의, 도킹, scoring, 감사추적)
- models/ : 독성 예측 baseline 모델
- tests/ : pytest 회귀 테스트
- docs/ : 한계 기록, 개발 진행상황, 발표 자료
- notebooks/ : Colab 노트북(번호순 개발 이력)

## 실행 환경

- Google Colab (GPU 런타임은 covalent docking 실험 시에만 필요)
- RDKit, ChEMBL API, AutoDock Vina, sascorer, Qwen(qwen3.8-max via
  DashScope) 등 외부 의존성은 각 노트북 세팅 셀에서 자동 설치

## 한계

docs/limitations.md에 상세 기록. 주요 항목: LLM 판단의 비결정성,
그리디 휴리스틱의 이론적 정당화 미흡, 일부 규칙(phosphor P=S 서브
클래스 등)은 치료용 약물 개선이라는 프로젝트 목적 범위 밖으로 판단해
제외.

## 라이선스 및 출처

RDKit, AutoDock Vina, sascorer(RDKit Contrib), ChEMBL/Guide to
Pharmacology 데이터를 활용. sascorer.py/fpscores.pkl.gz는 외부
저작물이라 저장소에 포함하지 않고 각 세션에서 다운로드.

## Contact

질문이나 피드백은 언제든 환영합니다.
- Email: hyekyeong.w@gmail.com
- GitHub: [@Dec32th](https://github.com/Dec32th)
