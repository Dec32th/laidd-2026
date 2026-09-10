# 미팅 심화 자료 — 구현 상세 포함

professor_meeting_prep.md의 요약을 뒷받침하는 실제 구현 상세.
교수님이 '어떻게 구현했는지' 물으실 때 바로 펼쳐서 보여줄 용도.

---

## 1. 전체 아키텍처

```
분자 입력
  -> detect_toxicophores()          [BRENK/PAINS 규칙 기반 진단]
  -> 다중문제 충돌조정                [남는 문제 수 최소화하는 순서로 정렬]
  -> propose_fix()                   [규칙별 치환 라이브러리에서 candidate 생성]
  -> 파괴적 편집 가드                 [원자손실 50%+ candidate 거부]
  -> should_debate()                 [rationale에 '[참고]' 있으면 토의 트리거]
  -> ask_llm_debate_fix()            [proposer-critic 다라운드 토의]
       ㄴ _try_get_docking_evidence  [등록된 표적이면 AutoDock Vina 자동 실행]
       ㄴ _try_compute_score         [QED/Lipinski/PAINS/SA/독성 종합점수]
       ㄴ _try_get_activity_risk     [Tanimoto/3D형태 활성보존 위험도]
       ㄴ _try_get_warhead_reference [문헌기반 워헤드 반응성 참고표]
  -> generate_audit_report()         [모든 판단과 근거를 사람이 읽을 리포트로]
```

## 2. 핵심 사례별 구현 상세

### 2.1 Aliphatic_long_chain 재설계

**문제**: BRENK의 실제 SMARTS(`[R0&D2][R0&D2][R0&D2][R0&D2]`)는 원소 종류와
무관하게 '고리 밖 degree-2 원자 4개 연속'이라는 순수 위상학적 패턴.
따라서 에테르 산소를 삽입해도(산소도 degree-2라서) 매치가 안 없어짐.

**해결**: `atom_editor.py`의 `insert_atom_multi_chain`에서 메틸 분기
(degree-2 -> degree-3)로 전략 전환. 양방향 체인 추적, 여러 매치 지점
탐색, `GetTotalNumHs()>0` 가드로 안전한 분기 위치만 선택.

**검증**: 규칙 미비 시 vs 수정 후 A/B 비교, Aliphatic_long_chain 차단
140건 -> 21건, 다른 규칙 회귀 없음, success 794 -> 903.

### 2.2 Michael_acceptor_1 EGFR 오귀속 편향

**문제**: `DOCKING_TARGETS`에 Michael_acceptor_1의 유일한 등록 표적이
EGFR(6JX4)이다 보니, 이 규칙이 걸리는 모든 분자가 무조건 EGFR에
도킹되고, critic이 도킹 점수만 보고 '이 분자는 EGFR 표적이다'라는
근거 없는 서사를 매번 자동 생성.

**해결**: `DOCKING_TARGETS` 각 항목에 `caveat` 필드 추가.
표적 특이성 근거가 없는 경우 명시적 경고 문구를 critic 프롬프트에
삽입(`agent.py`의 `_try_get_docking_evidence`가 caveat을 함께 반환).
캐시가 옛 caveat 텍스트를 그대로 반환하던 부수적 버그도 함께 수정
(점수는 캐시 재사용, caveat은 항상 최신 반영).

**검증**: 재현 5건 모두 critic이 'EGFR 벤치마크일 뿐'이라고 caveat을
인지하며 판단 근거를 도킹 수치에서 precedent 자체의 논리로 옮김.
200개 규모 재검증에서 판정이 3/3 반려 -> 3/3 승인으로 역전 확인.

### 2.3 공유결합 도킹 — 실패 두 번, 성공 한 번

**시도 1 (실패)**: RDKit SMILES 파싱 순서로 얻은 워헤드 원자 인덱스를
Vina 출력 PDBQT 좌표 리스트에 그대로 재사용. OpenBabel이 PDB->PDBQT
변환 시 원자 순서를 재배열해 완전히 다른 원자를 가리킴(원소 종류
일치율 20/37=54%, 우연 수준으로 확인).

**시도 2 (실패)**: 인덱스 매핑 대신 PDBQT를 mol로 재구성해 그 안에서
직접 워헤드 SMARTS 검색. PDBQT가 결합차수 정보를 온전히 담지 않아
OpenBabel이 C=C 이중결합을 올바르게 추측 못 함(실제 공유결합
억제제인 오시메르티닙으로도 재현 실패).

**성공 (Meeko 재도전)**: Meeko는 결합차수 손실 없는 PDBQT<->RDKit
왕복 변환을 공식 지원. Tethered Docking 방식(`--tether_smarts`,
`--rec_residue`)으로 EGFR Cys797(SG 좌표 실측: -52.745, -6.415,
-17.849, 6JX4 A chain)과 리간드를 물리적으로 연결한 상태로 도킹.

주요 구현 단계(scripts/meeko_covalent_docking_reproduce.py에 정리):
1. `mk_prepare_receptor.py`로 원본 리간드(YY3) 제거, 결측 곁사슬
   `--allow_bad_res`로 자동 처리, `-f`로 Cys797을 flexible로 분리
2. `mk_prepare_ligand.py`로 SDF 형태 리간드에 `--tether_smarts`로
   워헤드 지정, receptor의 Cys797과 연결된 tethered PDBQT 생성
3. `autogrid4`로 rigid receptor의 그리드 맵 계산
4. AutoDock-GPU(`--flexres`에 tethered 파일 전달)로 실제 도킹

**검증 결과**: 오시메르티닙 -5.31 kcal/mol(안정) vs 이타콘산
+7.06 kcal/mol(매우 불안정), 12.4 kcal/mol의 뚜렷한 차이로 방법의
변별력 확인.

**통합 보류 이유**: Colab 런타임이 매번 초기화되어 AutoDock-GPU를
세션마다 재컴파일해야 함(정적 바이너리가 아니라 소스 컴파일).
상시 파이프라인 통합의 운영비용이 크다고 판단해 보류.

### 2.4 LLM 비결정성

**발견**: 동일한 200개 샘플에 토의 파이프라인을 반복 실행했을 때
완전해결 수치가 67.2%~69.7%로 갈림 — 같은 입력에도 critic 판정이
매번 다를 수 있음을 확인.

**대응**: `ask_llm_debate_fix_consistent`(agent.py) — 같은 토의를
N회(기본 3회) 반복해 다수결로 최종 판정. `iterative_fix_loop`/
`batch_iterative_fix_loop`에 `use_consistent_debate` 파라미터로
연결(기본값 False로 기존 동작 보존).

**한계**: 50개 소규모 샘플로 효과를 검증하려 했으나, 토의 자체가
트리거된 횟수가 양쪽 다 3건뿐이라 통계적으로 유의미한 비교가
안 됨(표본 부족). 기능은 완성했지만 대규모 검증은 시간/비용 대비
효과가 불확실해 보류, 기본값 False 유지.

### 2.5 파괴적 편집 가드

**문제**: 다중문제 충돌조정이 '남는 문제 수'만 보다 보니, 분자를
통째로 파괴(원자 최대 93% 손실)해서 우연히 단순해진 candidate를
'최선'으로 착각하는 경우 발견. 예: 포레이트류 농약 분자에서
het-C-het_not_in_ring의 remove_substituent가 분자 대부분을 삭제,
C=S 하나만 남기고 success 처리.

**해결**: `molecule_editor.py`에서 candidate 채택 전 원자 손실 비율
계산, `destructive_edit_threshold`(기본 0.5) 이상이면 거부하고 다음
candidate로 재시도.

**민감도 분석**: 임계값을 0.3~0.7로 바꿔가며 300개 샘플 재실행,
5개 값 전부 동일한 결과(success 214/300). 원인: 실제 candidate
원자손실 비율 분포가 95.3% 30% 미만에 몰린 양극화 분포라 회색지대가
거의 없음 — '0.5는 임의값이지만 결과 민감도는 낮다'는 결론.

### 2.6 다중문제 충돌조정

**구현**: 여러 toxicophore가 동시에 있을 때, 각 규칙을 먼저 고쳤을
때 남는 전체 문제 수가 가장 적어지는 순서로 그리디 정렬
(`conflict_resolution='greedy'`, 기존 '리스트 순서대로'는
`'list_order'`로 실험 전환 가능하게 파라미터화).

**검증**: 300개 샘플에서 그리디 214 vs 리스트순서 213, 차이는 1건.
그 1건을 감사추적으로 직접 추적한 결과, 그리디가 '막다른 골목으로
이어지는 경로'를 정확히 회피해 성공한 교과서적 사례로 확인.
다중문제 동시발생 자체가 드물어 이 규모에서 개선폭은 작음(0.33%p).

## 3. 오픈소스 공개 관련 조언 요청 목록

- 지금 구조(RDKit 규칙 + LLM API 의존)로 오픈소스 공개 시, Qwen
  API 키가 없는 사용자를 위한 대체 경로(예: 규칙 기반만 동작)를
  어느 수준까지 지원해야 적절한지
- sascorer.py/fpscores.pkl.gz처럼 외부 저작물을 의도적으로 저장소에
  포함하지 않고 매 세션 다운로드하는 현재 방식이 라이선스 측면에서
  적절한지, 더 나은 관행이 있는지
- 연구 재현성 관점에서 노트북 기반 개발 이력(45~62번)을 그대로
  공개하는 것과, 정리된 스크립트/패키지 형태로 재구성하는 것 중
  어느 쪽이 대학원 지원 포트폴리오로 더 나은지
- Tox21처럼 널리 쓰이는 벤치마크 데이터셋 기반 도구를 공개할 때
  주의해야 할 점(예: 데이터 사용권, 인용 방식)
- 논문화를 염두에 둔다면, 지금 수준의 방법론적 검증(민감도분석,
  A/B 비교, 정직한 실패 기록)이 어느 정도까지 요구되는지

## 4. 이 모델을 한 문단으로 설명한다면

Tox21 데이터셋 분자의 독성 구조를 BRENK/PAINS 규칙으로 진단하고,
화학 지식 기반 치환 라이브러리에서 candidate를 생성한 뒤, 그 중
실제 승인약물 사례와 유사할 가능성이 있는 candidate에 한해 LLM
proposer-critic 토의를 거치게 하는 시스템입니다. 이 토의는 단백질
도킹 점수, 합성용이성, 활성보존 지표, 문헌 기반 워헤드 반응성
정보를 근거로 받아 판단하며, 모든 판단 과정(승인/반려/사람에게
위임)을 감사추적 형태로 기록합니다. 핵심은 'LLM이 화학구조를 직접
생성하지 않고, 결정론적 화학 도구가 만든 후보들 중에서만 판단한다'는
설계 원칙이며, 개발 과정 자체가 여러 차례의 가설-검증-반증-재도전을
거친 반복적 탐구였습니다.