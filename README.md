# AI 기반 신약 후보물질 독성 최적화 에이전트

제4회 AI 신약개발 경진대회(JUMP AI 2026) 참가 프로젝트이며, 동시에
2026년 한 해 동안 진행하는 LAIDD 학습 기록을 담은 리포지토리입니다.

## 무엇을 하는 프로젝트인가

RDKit과 Tox21/Ames/hERG/DILI 데이터를 활용해, 독성이 우려되는 화합물을
**예측 → 구조 진단 → 치환 판단(LLM) → 구조 개선 → 활성 보존 검토 →
재평가**하는 과정을 자동으로 반복하는 에이전틱 최적화 루프를 만듭니다.

1. Tox21/Ames/hERG/DILI 기반 모델로 화합물의 독성 예측
2. RDKit FilterCatalog로 문제 구조(toxicophore)를 탐지
3. LLM 에이전트가 (여러 문제 중) 무엇을 먼저 고칠지, (여러 후보 중) 어떤
   치환이 화학적으로 타당한지 판단 — 단, 실제 구조 변경은 RDKit 기반
   결정적 도구가 전담(LLM은 화학구조를 직접 생성하지 않음)
4. 3개 전문가 에이전트(독성학/의약화학/약리학)가 정량 지표(2D/3D 구조
   유사도, QED, LogP, 합성용이성)를 근거로 활성 보존 가능성을 함께
   검토, 불확실한 경우 사람 검토를 명시적으로 요청
5. 개선된 구조를 다시 예측 모델에 넣어 독성이 실제로 줄었는지 확인,
   문제가 남으면 반복

## 현재 상태 (2026-08 기준)

- 치환 라이브러리 35개 규칙, 11가지 이상의 구조 편집 방식
- Held-out(valid set) 커버리지 33.2%, 단일 문제 분자 기준 최소 1단계
  이상 개선 84%
- 4개 독립 독성 endpoint(Tox21/Ames/hERG/DILI) 교차검증 (Ames/Tox21/
  DILI 유의, hERG 미유의)
- 3-에이전트(독성학/의약화학/약리학) 판단 + 활성보존 근사지표 통합,
  실제 승인·철수 약물(ChEMBL 23건) 기반 소급 검증
- COMT 표적 도킹 검증(스트레치 목표): catechol 규칙의 rationale을
  실제 결합 스코어로 확인(도파민 -5.72 vs 메톡시-도파민 -5.41 kcal/mol)
- 개발 과정에서 발견한 이름 불일치·중복 규칙·재시도 로직 결함 등의
  오류를 전수검사 및 ablation 실험으로 스스로 발견·수정 (자세한 내용은
  `docs/limitations.md` 참고)

## 폴더 구조

- `src/tools/` — 핵심 로직 (독성 진단, 치환 라이브러리, 구조 편집,
  LLM 에이전트 판단, 시각화)
- `notebooks/` — 개발 과정 노트북 (번호 순, 각 노트북 상단에 진행상황 요약)
- `models/` — 학습된 baseline 모델 메타데이터
- `outputs/` — 실험 결과 로그(JSON/텍스트), 재현 가능한 산출물
- `docs/` — 한계 기록, 설계 결정, 제안서 첨부용 검증 문서(ablation,
  case study, 합성가능성, 지표-LLM 상호보완성, 실험결과 로그, 로드맵 등)

## 무엇을 할 것인가

- **공모전 준비 (~8/7)**: 제안서 최종 정리, hwpx 포맷팅, test set 최종
  검증
- **2026년 전반**: LAIDD 강의를 학기 일정과 병행해 꾸준히 수강하며, 이
  리포지토리에 학습 내용을 정리
- **공모전 이후**: 선례 라이브러리 확장, 추가 표적 도킹 검증, 특정
  치료 영역 특화 등 장기 포트폴리오로 확장 예정 (`docs/toxguard_roadmap.md`
  참고)

## 기술 스택

Python, RDKit, Tox21, Ames(TDC), hERG(TDC), DILI(TDC), MMPDB,
AutoDock Vina, OpenBabel, Google Colab, Gemini/Qwen(LLM 판단),
py3Dmol(3D 시각화)

## 참고

- [JUMP AI 공모전 공고](https://www.laidd.org/mod/ubboard/article.php?id=1&bwid=879)
- [LAIDD](https://www.laidd.org)
