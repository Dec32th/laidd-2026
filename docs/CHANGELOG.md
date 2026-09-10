# CHANGELOG

노트북 번호 순 개발 이력 요약. 예선 제출 이후에도 지속 개발 중.

## 예선 제출까지 (초기 개발)

- 규칙 라이브러리 33.2% -> 48.9% 커버리지 확장 (빈도 기반 규칙 추가:
  Aliphatic_long_chain, isolated_alkene, quaternary_nitrogen_1/2,
  phenol_ester 등)
- diketo_group 규칙 추가, phosphor/halogenated_ring_1/heavy_metal/
  Oxygen-nitrogen_single_bond/iodine/Perfluorinated_chain은 정당한
  근거로 보류 결정(과도하게 광범위하거나 승인약물에 흔하거나 산업
  화학물질 특이적)
- 3-agent 아키텍처, 도킹 검증 스텁, precedent_library.py 최초 도입
  (11건), 다중목적 scoring 설계

## 노트북 45-47: 버그 수정 및 선례 확장

- Aliphatic_long_chain SMARTS 실측 분석 후 재설계(에테르 삽입 ->
  분기 전략)로 완전해결
- PRECEDENT_LIBRARY 11건 -> 23건 확장 (ChEMBL/GtoPdb API 실조회
  기반, 각 항목 승인약물/기전 검증)

## 노트북 48-49: LLM 오케스트레이션 + 도킹 자동화

- LLM 호출 병렬화(batch_iterative_fix_loop), 토의 트리거 조건/예산
  상한 도입
- proposer-critic 다라운드 토의 로직(ask_llm_debate_fix) 구축
- docking.py: AutoDock Vina 자동화, 4개 표적(COMT, EGFR, NQO1 x2)
  검증. 도킹 근거를 critic 프롬프트에 연결

## 노트북 50-51: 다중문제 처리 + 안전장치

- 다중문제 충돌조정(그리디, 남는 문제 수 최소화) 도입
- 감사추적 리포트(audit.py) 구축, 이후 반려/에스컬레이트 케이스
  기록 공백 발견 및 수정
- pytest 회귀 테스트 스위트 도입
- Tox21 baseline 모델(tox_baseline.py), 활성보존 지표(activity_metrics.py)
  통합

## 노트북 52-53: 최종 측정 + 안정화

- 전체 valid set(1173개) 첫 공식 측정: 규칙기반 74.5%, 토의파이프라인
  67.2-69.7%(실행 간 변동 발견)
- LLM 비결정성 완화(ask_llm_debate_fix_consistent, 다수결)
- Michael_acceptor_1 도킹 caveat 강화, 캐시가 옛 caveat을 반환하던
  버그 수정

## 노트북 54: 공유결합 근사 시도와 정직한 실패

- 표준 도킹 결과에 원자 인덱스를 매핑해 공유결합 형성 가능성을
  근사하려 두 차례 시도 (SMILES-PDBQT 인덱스 매핑, PDBQT->mol
  재구성) - 둘 다 원자순서/결합차수 정보 손실로 실패, 정직하게 폐기
- 대안으로 워헤드 반응성 참고표(WARHEAD_REACTIVITY_REFERENCE) 도입,
  이타콘산의 실제 KEAP1/GAPDH 공유결합 문헌 사례 반영

## 노트북 55-57: 확장 및 재검증

- 워헤드 참고표를 alkyl_halide, disulphide로 확장
- phosphor 규칙의 P=S 서브클래스는 유기인계 농약 골격으로 확인,
  치료용 약물 개선 목적 범위 밖으로 null-result 기록
- 파괴적 편집 가드 임계값(0.3-0.7) 및 충돌조정 방식(그리디 vs
  리스트순서) 민감도 분석
- 전체 valid set 재측정으로 수치 안정성 확인

## 노트북 58: Meeko + AutoDock-GPU 기반 covalent docking 성공

- Meeko의 결합차수 보존 기능으로 노트북 54의 원자매핑 문제를
  근본적으로 우회
- Tethered docking 방식으로 EGFR Cys797-오시메르티닙 공유결합 형태의
  결합 에너지 계산 성공: 오시메르티닙 -5.31 kcal/mol(안정) vs
  이타콘산 +7.06 kcal/mol(불안정), 12.4 kcal/mol 뚜렷한 차이로
  방법 검증
- Colab 매 세션 재컴파일 비용 때문에 상시 파이프라인 통합은 보류,
  성공 사례로 기록

## 노트북 59-61: 마무리 및 문서화

- 워헤드 참고표 확장 대상 최종 확인 (epoxide/acid_halide는 해당 없음)
- Meeko 재현 스크립트 정리, 미팅용 감사추적 예시 3건 작성
- no_known_fix 재조사: iodine(덱스트로티록신/아미오다론 근거),
  het_thio_666_A(13)/페노티아진(9개 승인 항정신병약 근거) 선례 추가
- iodine의 예전 "100% 중복" 결론이 현재 데이터로는 성립하지 않음을
  재검증으로 확인
- README.md, CHANGELOG.md 작성

## 프로젝트 상태

JUMP AI 2026 공모전(4th JUMP AI competition, 팀명 MorForge) 예선
탈락 확정. 예선 결과와 무관하게 오픈소스 공개 및 대학원 진학
포트폴리오로 지속 개발 중.
