

"""선례 라이브러리 — 승인/철수 약물, 정량 활성 데이터, 도킹 검증 결과를
판단 에이전트 프롬프트에 실시간 주입하기 위한 구조화된 근거 저장소.
모든 항목은 이 세션에서 ChEMBL/GtoPdb API 조회 또는 실제 도킹 실행으로
직접 확인한 것만 포함한다(추정/일반 지식은 배제).
"""

PRECEDENT_LIBRARY = [
    {"rule": "Thiocarbonyl_group", "type": "긍정_승인약물쌍",
     "description": "티오펜탈(C=S)/펜토바비탈(C=O), 티아밀랄(C=S)/세코바비탈(C=O) - "
                     "동일 사이드체인, C=S->C=O만 다른 실제 승인 마취제 쌍. "
                     "baseline 모델 기준 옥소형이 티오형보다 Tox21 평균 예측값 낮음(-0.008~-0.009)."},
    {"rule": "catechol", "type": "정량_활성데이터",
     "description": "도파민이 D1(Ki 4.3-5.6nM)/D2(Ki 4.7-7.2nM)/D3(Ki 6.4-7.3nM) 수용체에 "
                     "단자릿수 nM 강력 결합 - 카테콜 골격이 활성에 필수적임을 정량적으로 뒷받침."},
    {"rule": "hydroxamic_acid", "type": "부정_참고사례_검증필요",
     "description": "하이드록삼산 골격(보리노스타트 등 HDAC 억제제)은 아연 킬레이션이 "
                     "약효 핵심이므로, 이 계열에 대한 무분별한 치환은 약효 상실 위험. "
                     "(문헌 재확인 필요)"},
    {"rule": "beta-keto/anhydride", "type": "긍정_통계검증결과",
     "description": "MMPDB 공식 통계 도구로 재검증한 결과, Tox21 규모(1173개)에서 "
                     "무수물 관련 매칭쌍은 표본 부족(count=1)으로 통계적 유의성 확보 불가 - "
                     "데이터형 접근보다 문헌형 근거가 더 신뢰할 만함을 시사."},
    {"rule": "Michael_acceptor_1", "type": "위험=메커니즘_참고",
     "description": "에타크린산(이뇨제, FDA 승인)은 시스테인 잔기와의 공유결합 자체가 "
                     "작용 메커니즘인 공유결합 억제제 - Michael acceptor 경고가 항상 "
                     "제거 대상은 아님을 보여주는 실제 승인약물 사례."},
    {"rule": "alkyl_halide", "type": "위험=메커니즘_참고",
     "description": "메클로르에타민, 사이클로포스파미드 등 알킬화 항암제는 DNA 알킬화 "
                     "반응성 자체가 세포독성 치료 메커니즘 - 이 계열에는 할로겐 제거가 "
                     "부적절함을 보여주는 실제 승인약물 사례."},
    {"rule": "azo_A(324)", "type": "위험=메커니즘_참고_검증완료",
     "description": "설파살라진(SMILES 내 /N=N/ 아조 결합 확인, ChEMBL max_phase=4.0, "
                     "GtoPdb FDA 승인 1950년/WHO 필수의약품)은 아조 결합이 장내 "
                     "세균에 의해 환원되어 활성 대사물(5-ASA)을 방출하는 프로드러그 - "
                     "실제 조회로 검증됨."},
    {"rule": "catechol", "type": "도킹검증_결과",
     "description": "COMT(PDB 1VID) 도킹 검증: 도파민(-5.72 kcal/mol)→메톡시도파민"
                     "(-5.41 kcal/mol), 변화폭 +0.31 kcal/mol로 약화 방향이나 이는 "
                     "1 kcal/mol 미만의 작은 차이로 도킹 자체의 오차범위 내일 수 있어 "
                     "단정적 근거로 삼기엔 약함. 에피네프린은 반대로 미세 강화"
                     "(-6.21→-6.32, -0.10) - 두 경우 모두 변화폭이 작아, 도킹 수치보다는 "
                     "카테콜의 수용체 결합 필수성(정성적 근거)이 더 강한 판단 기준."},
    {"rule": "Michael_acceptor_1", "type": "도킹검증_방법론한계",
     "description": "EGFR(PDB 6JX4) 도킹 검증: 오시메르티닙(-7.13)→C=C환원버전(-7.08), "
                     "거의 무변화(+0.05). 표준(비공유) 도킹이 오시메르티닙의 실제 "
                     "공유결합(Cys797) 메커니즘을 포착하지 못하는 방법론적 한계 확인 - "
                     "공유결합 억제제 계열은 일반 도킹 스코어만으로 활성 손실을 판단하지 "
                     "말 것(도킹 무변화가 곧 활성 유지를 뜻하지 않음)."},
    {"rule": "hydroquinone", "type": "도킹검증_결과",
     "description": "NQO1 도킹 검증: 퀴논(-3.29)→하이드로퀴논(-4.08), 결합 강화(-0.79, "
                     "1 kcal/mol에 근접하는 뚜렷한 변화). 메틸퀴논(-3.80)→환원버전"
                     "(-4.29)도 강화(-0.49), 2건 모두 일관되게 강화 방향. NQO1이 실제로 "
                     "퀴논을 하이드로퀴논으로 환원하는 효소이므로, 이 치환 방향은 해독 "
                     "반응경로와 자연스럽게 정렬되며 실측 결합력도 개선됨 - 활성 손실 "
                     "우려가 낮은 것으로 확인됨."},
    {"rule": "quinone_A(370)", "type": "도킹검증_결과",
     "description": "NQO1 도킹 검증: 퀴논(-3.29)→하이드로퀴논(-4.08), 결합 강화(-0.79). "
                     "실제 표적 효소와의 결합력이 오히려 개선되는 것으로 실측 확인됨 "
                     "(hydroquinone 규칙과 동일 표적 데이터 공유)."},
    {"rule": "Aliphatic_long_chain", "type": "긍정_승인약물_확인(구조는_동의어로_대체확인)",
     "description": "POLIDOCANOL(라우릴알코올+에틸렌옥사이드 평균 9개 반복부가체)은 ChEMBL 조회로 "
                     "승인 확인됨(max_phase=4.0, first_approval=2010, ATC C05BB02, "
                     "dosed_ingredient=True, withdrawn=False, 상품명 Asclera/Aethoxysklerol). "
                     "ChEMBL에 단일 SMILES는 없으나(polymer_flag=1, structure_type=NONE) "
                     "이는 다분산 고분자라 원천적으로 단일 구조가 없기 때문이며, 공식 동의어"
                     "(USP: Polyoxyl 9 lauryl ether, JAN: Lauromacrogol 400)가 "
                     "\"장쇄 알킬+반복 에테르\" 구조를 명확히 정의함 - 실제 승인약물에서 "
                     "이 전략이 쓰이고 있음을 뒷받침."},
    {"rule": "Aliphatic_long_chain", "type": "부정_참고사례_검증필요",
     "description": "ChEMBL 서브구조 검색(에테르 삽입 사슬 모티프)으로 매치된 승인약물은 "
                     "에리스로마이신/아지스로마이신/암포테리신B였으나, 매치 위치를 IsInRing으로 "
                     "확인한 결과 전부 매크로락톤/당 고리 내부의 고리형 에테르로, 우리 규칙이 "
                     "다루는 \"고리 밖 열린 사슬\" 상황과는 구조적으로 다름 - 이 계열은 직접적 "
                     "근거로 부적합함이 확인됨."},
    {"rule": "isolated_alkene", "type": "긍정_승인약물쌍",
     "description": "SIROLIMUS(시롤리무스), TACROLIMUS ANHYDROUS(타크로리무스) - 둘 다 ChEMBL "
                     "조회로 승인·비철수 확인됨(withdrawn_flag=False), 대형 매크로라이드 면역억제제로 "
                     "현재도 널리 처방됨. problem_smarts로 직접 매치되는 고립 지방족 알켄이 구조 "
                     "안에 실제 존재 - 고립 알켄이 항상 제거 대상은 아님을 보여주는 실제 승인약물 사례."},
    {"rule": "isolated_alkene", "type": "위험=메커니즘_참고_인과불명",
     "description": "CYCLOBARBITAL, HEXOBARBITAL 둘 다 ChEMBL 조회로 withdrawn_flag=True 확인됨, "
                     "둘 다 problem_smarts에 매치되는 사이클로헥세닐 고립 알켄 치환기를 가짐. "
                     "다만 바르비투르산염 계열은 호흡억제·의존성 등 일반적 안전성 문제로 철수된 "
                     "사례가 많아, 이 알켄 구조가 철수의 직접 원인이라는 인과관계는 확인되지 않음 "
                     "(상관관계만 관찰, 문헌 추가 확인 필요)."},
    {"rule": "nitro_group", "type": "위험=메커니즘_참고_검증완료",
     "description": "METRONIDAZOLE, NITROFURANTOIN, BENZNIDAZOLE 셋 다 ChEMBL 조회로 승인·비철수 "
                     "확인됨(max_phase=4.0, withdrawn_flag=False), SMILES에 니트로기([N+](=O)[O-]) "
                     "실제 존재 확인. 항균/항기생충제 계열에서 니트로기의 선택적 환원 활성화 자체가 "
                     "치료 메커니즘인 프로드러그 설계 사례 - 이런 계열에는 니트로기 제거가 "
                     "부적절함을 실제 조회로 검증함."},
    {"rule": "aniline", "type": "위험=메커니즘_참고_검증완료",
     "description": "SULFANILAMIDE, SULFAMETHOXAZOLE, PROCAINAMIDE 셋 다 ChEMBL 조회로 승인·비철수 "
                     "확인됨(max_phase=4.0, withdrawn_flag=False), SMILES 확인 결과 셋 다 아실화되지 "
                     "않은 유리 1차 방향족 아민(아닐린) 형태로 실제 처방됨. 설파계 항생제·항부정맥제 "
                     "계열에서 특이체질 반응 위험에도 불구하고 유리 아닐린 골격이 오랜 기간 널리 "
                     "쓰여온 사례 - 이 경고가 절대적 배제 기준이 아님을 실제 조회로 검증함."},
    {"rule": "Sulfonic_acid_2", "type": "위험=메커니즘_참고_검증완료",
     "description": "LISDEXAMFETAMINE DIMESYLATE, SAQUINAVIR MESYLATE 둘 다 ChEMBL 조회로 승인·비철수 "
                     "확인됨(max_phase=4.0, withdrawn_flag=False). RDKit GetMolFrags로 분자 조각을 "
                     "분리해 확인한 결과, 설폰산(메실산) 매치는 둘 다 작은 카운터이온 조각(CS(=O)(=O)O, "
                     "5원자)에서만 나오고 주 약효 골격 조각(19원자, 49원자)에서는 전혀 매치되지 않음 - "
                     "설폰산이 활성 골격이 아니라 순수 염 형성용 카운터이온인 경우가 실제로 존재함을 "
                     "구조적으로 검증함. 이런 경우 본 규칙의 치환 대상이 아님."},
    {"rule": "phosphor", "type": "위험=메커니즘_참고_검증완료",
     "description": "FOSPHENYTOIN(유리산 형태) ChEMBL 조회로 승인·비철수 확인됨(max_phase=4.0, "
                     "withdrawn_flag=False), problem_smarts 실제 매치 확인. 페니토인의 인산에스터 "
                     "프로드러그로, 체내 인산가수분해효소에 의한 에스터 절단 자체가 설계된 방출 "
                     "메커니즘 - 이 경우 인산에스터 절단이 규칙이 우려하는 신경독성 반응성이 아니라 "
                     "오히려 활성화 경로이므로, 프로드러그 맥락에서는 본 규칙의 무분별한 적용이 "
                     "부적절할 수 있음을 실제 조회로 검증함."},
]


def get_precedents(rule_name: str) -> str | None:
    """규칙 이름으로 관련 선례를 찾아 프롬프트에 넣을 텍스트로 반환."""
    matches = [p for p in PRECEDENT_LIBRARY if p['rule'] == rule_name]
    if not matches:
        return None
    return "\n".join([f"- [{m['type']}] {m['description']}" for m in matches])
