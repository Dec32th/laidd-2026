
REPLACEMENT_LIBRARY = {
    "nitro_group": {
        "problem_smarts": "[N+](=O)[O-]",
        "candidates": [
            {"smiles": "N", "name": "primary amine",
             "rationale": "[참고] 메트로니다졸, 니트로푸란토인, 벤즈니다졸 등 일부 "
                          "항균제/항기생충제는 니트로기의 선택적 환원 활성화 자체가 "
                          "치료 메커니즘이므로, 이런 프로드러그 설계 맥락에서는 본 "
                          "치환이 적절하지 않을 수 있음. || 극성을 유지하면서 니트로기의 "
                          "환원성 대사 중간체 생성 경로를 제거함"},
            {"smiles": "S(=O)(=O)N", "name": "sulfonamide",
             "rationale": "약물유사 골격에서 흔히 쓰이는 안정적 대체기로, 수소결합 donor/acceptor 특성을 일부 유지"},
            {"smiles": "C#N", "name": "nitrile",
             "rationale": "대사 안정성이 개선된 사례가 문헌에 다수 보고됨, 다만 극성은 다소 감소"},
        ],
    },
    "aldehyde": {
        "problem_smarts": "[CX3H1](=O)",
        "candidates": [
            {"smiles": "C(=O)N", "name": "amide",
             "rationale": "알데히드의 친전자성(단백질 부가물 형성 우려)을 제거하면서 유사한 형태 유지"},
            {"smiles": "C(O)", "name": "alcohol",
             "rationale": "가장 단순한 환원형 대체, 반응성 크게 감소"},
        ],
    },
    "Michael_acceptor_1": {
        "edit_method": "atom_edit",
        "problem_smarts": "C=CC(=O)",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "saturated (C-C single bond)",
             "rationale": "[참고] 에타크린산처럼 시스테인 잔기와의 공유결합 자체가 "
                          "작용 메커니즘인 공유결합 억제제(covalent inhibitor) "
                          "계열에는 본 경고가 그대로 적용되지 않을 수 있음. || "
                          "알파,베타-불포화 카르보닐의 C=C 이중결합을 환원하여 "
                          "단백질 친전자성 부가반응(Michael addition, covalent "
                          "binding) 위험을 제거함"},
        ],
    },
    "acid_halide": {
        "problem_smarts": "C(=O)[F,Cl,Br,I]",
        "candidates": [
            {"smiles": "C(=O)N", "name": "amide",
             "rationale": "고반응성 아실할라이드를 안정적인 아마이드로 대체"},
            {"smiles": "C(=O)O", "name": "ester",
             "rationale": "아마이드보다 극성이 낮고 유연한 대체 옵션, 가수분해 속도 조절 가능 (검증 필요)"},
        ],
    },
    "alkyl_halide": {
        "problem_smarts": "[Cl,Br,I]",
        "candidates": [
            {"smiles": "O", "name": "hydroxyl (alcohol)",
             "rationale": "[참고] 메클로르에타민, 사이클로포스파미드, 카머스틴, "
                          "클로람부실 등 알킬화 항암제는 DNA 알킬화(반응성) 자체가 "
                          "세포독성 치료 메커니즘이므로, 이 계열에는 본 치환이 "
                          "적절하지 않음. || 이탈기를 제거해 알킬화 반응성을 없앰, "
                          "극성은 유사하게 유지"},
            {"smiles": "F", "name": "fluorine",
             "rationale": "할로겐을 유지하되 C-F 결합은 강해 이탈기로 작용하지 않음, 입체적 크기도 유사"},
        ],
    },
    "aniline": {
        "edit_method": "atom_edit",
        "problem_smarts": "[NH2]c1ccc([#6,#7,#8,#16])cc1",
        "target_idx_in_pattern": 0,
        "ring_atom_indices_in_pattern": [1, 2, 3, 4, 6, 7],
        "anchor_indices_in_pattern": (0, 5),
        "candidates": [
            {"edit_type": "add_substituent", "param": "C(=O)C",
             "target_idx_in_pattern": 0,
             "name": "acetamide (acylated amine)",
             "rationale": "[참고] 설파계 항생제(설파닐아마이드, 설파메톡사졸 등)와 "
                          "프로카인아마이드처럼 아닐린 골격이 반응성 대사가 아닌 "
                          "안정적 형태로 널리 처방되어 온 사례가 다수 있음. 이 경우 "
                          "특이체질 반응은 드물고 예측이 어려워, 본 경고를 절대적 "
                          "배제 기준이 아닌 참고 신호로 해석해야 함. || 1차 방향족 "
                          "아민을 아마이드로 아실화하여 N-hydroxylation 경로 자체를 차단"},
            {"edit_type": "replace_ring", "param": "[*:1]C12CC(C1)(C2)[*:2]",
             "ring_atom_indices_in_pattern": [1, 2, 3, 4, 6, 7],
             "anchor_indices_in_pattern": (0, 5),
             "name": "BCP (bicyclo[1.1.1]pentane)",
             "rationale": "para-이치환 아닐린의 방향족 벤젠 고리를 포화 bicyclic "
                          "탄소골격(BCP)으로 교체함. 방향족성 제거로 aniline reactive "
                          "metabolite(RM) 형성 및 CYP-inhibition을 감소시켜, 퀴논이민 "
                          "생성 경로를 차단하고 특이체질 약물 부작용(IADR) 위험을 낮춤 "
                          "(문헌 근거, 학생 제공). 벤젠과의 공간적 유사성, Fsp3 증가, "
                          "실제 성공 사례가 많아 채택. 아마이드화(단순 아민 치환)보다 "
                          "변화 폭이 크지만, 물성 개선 효과도 더 큼"},
        ],
    },
    "Sulfonic_acid_2": {
        "problem_smarts": "S(=O)(=O)[OX2H1,OX1-]",
        "candidates": [
            {"smiles": "S(=O)(=O)N", "name": "sulfonamide",
             "rationale": "[참고] 암페타민 설페이트, 사퀴나비르 메실레이트처럼 "
                          "일부 승인약물에서 설폰산/설폰산 유사기는 활성 골격이 "
                          "아니라 염(salt) 형성을 위한 카운터이온으로만 존재함. "
                          "이 경우 본 규칙이 다루는 '독성 유발 골격'과 무관하므로, "
                          "치환 대상 여부를 판단하기 전에 이 산이 활성 골격의 "
                          "일부인지 염 형성용인지 구분이 필요함. || 생리적 pH에서 "
                          "이온화 정도(전하)를 크게 낮춰 세포막 투과성을 개선함. "
                          "설폰산은 대부분 음이온 상태로 존재해 경구 흡수가 저해되는 "
                          "경우가 많으나, 설폰아마이드는 유사한 골격을 유지하면서도 "
                          "중성에 가까워 약물유사성이 개선됨"},
            {"smiles": "C(=O)O", "name": "carboxylic acid",
             "rationale": "설폰산보다 산성도가 약하고 부피가 작은 산성 bioisostere "
                          "(검증 필요)"},
        ],
    },
    "imine_1_oxime": {
        "edit_method": "atom_edit",
        "problem_smarts": "C=N[OX2H1]",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "amine (reduced)",
             "rationale": "옥심의 C=N 결합을 환원하여, 가수분해 시 원래의 반응성 "
                          "카르보닐(알데히드/케톤)로 되돌아갈 수 있는 대사 불안정 "
                          "경로를 제거함"},
        ],
    },
    "imine_1_general": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CX3;!$(C(N)(N)=N)]=N",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "amine (reduced)",
             "rationale": "일반 이민(C=N-R)을 환원하여 가수분해 시 반응성 카르보닐로 "
                          "되돌아갈 수 있는 대사 불안정 경로를 제거함. 옥심 특유의 "
                          "메커니즘보다는 근거가 다소 약하며, 하위 구조별 개별 검증 필요. "
                          "구아니딘(N-C(=N)-N, 공명구조로 일반 이민과 반응성이 다름)은 "
                          "이 SMARTS에서 명시적으로 제외함"},
        ],
    },
    "catechol": {
        "edit_method": "atom_edit",
        "problem_smarts": "[OX2H;$(Oc1ccccc1O)]",
        "target_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "add_substituent", "param": "C", "name": "methoxy",
             "rationale": "[참고] 도파민, 에피네프린, 이소프로테레놀 등 카테콜아민류 "
                          "약물은 카테콜 구조 자체가 아드레날린/도파민 수용체 결합에 "
                          "필수적인 약효 골격이므로, 이 경우 본 치환은 독성 감소가 "
                          "아니라 약효 상실로 이어짐. || 인체의 COMT(catechol-O-"
                          "methyltransferase) 효소가 카테콜을 메톡시페놀로 메틸화하여 "
                          "해독하는 생리적 경로와 동일한 원리. 오르토-퀴논으로의 산화 "
                          "경로를 차단하여 세포독성/유전독성 우려를 낮춤 (학생 확인 "
                          "예정: ScienceDirect catechol overview, PMC6643002 등 참고)"},
        ],
    },
    "Thiocarbonyl_group": {
        "edit_method": "atom_edit",
        "problem_smarts": "[#6]=[#16]",
        "target_idx_in_pattern": 1,
        "candidates": [
            {"edit_type": "replace_element", "param": 8, "name": "carbonyl (O replacing S)",
             "rationale": "[참고] 티오펜탈·티아밀랄(치오바르비투레이트, C=S가 지용성 "
                          "증가로 빠른 마취효과에 기여)과 티오구아닌(퓨린 유사 항대사물, "
                          "황이 작용기전에 필수)처럼 황 원자가 약효/효력에 직접 "
                          "기여하는 경우가 있어, 이 계열에는 본 치환이 부적절할 수 "
                          "있음. || 황을 산소로 대체(티오카르보닐->카르보닐)하는 것은 "
                          "흔한 bioisostere 전략으로, 갑상선 기능 저해 등 황 함유 "
                          "작용기 특유의 대사/독성 우려를 낮춤 (검증 필요, "
                          "thiourea->urea 치환 논리와 동일 계열)"},
        ],
    },
    "thiol_2": {
        "problem_smarts": "[SX2H1]",
        "candidates": [
            {"smiles": "O", "name": "hydroxyl (alcohol)",
             "rationale": "티올의 금속 킬레이팅 및 산화(이황화물/술펜산 형성) 반응성을 "
                          "제거하면서, 극성·수소결합 특성을 유사하게 유지함"},
            {"smiles": "C(=O)N", "name": "amide",
             "rationale": "티올을 아마이드로 대체하여 반응성을 낮추면서 약물유사 골격에서 "
                          "흔히 쓰이는 안정적 작용기로 전환 (검증 필요)"},
        ],
    },
    "thiol_1": {
        "edit_method": "atom_edit",
        "problem_smarts": "C(=S)[SX1-]",
        "candidates": [
            {"edit_type": "replace_multi",
             "param": [
                 {"idx_in_pattern": 1, "new_element": 8, "new_charge": 0},
                 {"idx_in_pattern": 2, "new_element": 7, "new_charge": 0},
             ],
             "name": "carbamate (O,N replacing S,S)",
             "rationale": "디티오카바메이트(R-O-C(=S)-S-)를 카바메이트(R-O-C(=O)-N)로 "
                          "전환. 두 황 원자를 각각 산소·질소로 교체하여 금속 킬레이팅 "
                          "능력과 효소 억제 활성(디티오카바메이트류 특유의 살충제성 "
                          "독성 기전)을 제거함 (검증 필요)"},
        ],
    },
    "het-C-het_not_in_ring": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CX4](O)(O)",
        "candidates": [
            {"edit_type": "remove_substituent",
             "center_idx_in_pattern": 0,
             "remove_idx_in_pattern": 1,
             "upgrade_bond_to_idx_in_pattern": 2,
             "name": "ketone/ester (one alkoxy removed, C=O formed)",
             "rationale": "아세탈/케탈 또는 오르토에스터(탄소 하나에 알콕시기 2개 "
                          "이상)는 가수분해에 민감하여 반응성 카르보닐(케톤/알데히드)로 "
                          "쉽게 분해되며 대사 불안정성을 일으킴. 알콕시기 하나를 제거하고 "
                          "남은 산소를 카르보닐로 승격시켜, 가수분해로 어차피 도달할 "
                          "안정한 최종 형태로 미리 전환함 (검증 필요)"},
        ],
    },
    "hydroquinone": {
        "edit_method": "atom_edit",
        "problem_smarts": "[OX2H]c1ccc([OX2H,NX3H1,NX3H2])cc1",
        "target_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "add_substituent", "param": "C", "name": "methoxy",
             "rationale": "[참고] 아세트아미노펜은 정상 용량에서는 안전하며 과다복용 "
                          "시에만 위험한 용량 의존적 사례임. 본 시스템은 치료지수를 "
                          "고려하지 않으므로, 아트로핀·디곡신·와파린처럼 좁은 치료지수를 "
                          "가진 기존 약물 전반에 유사하게 적용되는 한계임. || 파라 "
                          "위치에 OH와 (OH 또는 NH)가 있는 구조(하이드로퀴논/파라-"
                          "아미노페놀 계열)는 산화되어 파라-퀴논 또는 파라-퀴논이민(예: "
                          "아세트아미노펜의 NAPQI)을 형성, 글루타치온 고갈과 단백질 "
                          "공유결합을 통한 간독성 위험이 있음"},
        ],
    },
    "azo_A(324)": {
        "edit_method": "atom_edit",
        "problem_smarts": "N=N",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "hydrazine (reduced)",
             "rationale": "아조기(N=N)는 체내에서 아조환원효소에 의해 환원되어 두 개의 "
                          "방향족 아민으로 분해되며, 그 중 일부(벤지딘류 등)가 발암성을 "
                          "가지는 것으로 잘 알려짐(아조 색소의 대표적 독성 메커니즘). "
                          "이중결합을 환원하여 하이드라진 형태로 전환, 완전한 아민 "
                          "분해 경로 자체를 차단함 (검증 필요: 하이드라진 자체의 "
                          "잔여 반응성은 추가 확인 필요)"},
        ],
    },
    "Three-membered_heterocycle": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CX4]1[OX2][CX4]1",
        "candidates": [
            {"edit_type": "open_epoxide", "break_pair_in_pattern": (1, 2),
             "name": "vicinal diol (ring-opened)",
             "rationale": "에폭시드(3원자 고리, 옥시란)는 고리 변형(strain)으로 인해 "
                          "친핵체(DNA, 단백질)와 쉽게 반응하는 알킬화제로 작용함. "
                          "체내 에폭시드 가수분해효소(epoxide hydrolase)가 실제로 "
                          "수행하는 반응과 동일하게 고리를 열어 비시날 디올(vicinal "
                          "diol)로 전환, 반응성을 제거함"},
        ],
    },
}

def get_replacement_candidates(rule_name: str) -> dict | None:
    """rule_name에 해당하는 치환 정보(SMARTS + 후보 리스트)를 반환. 없으면 None."""
    return REPLACEMENT_LIBRARY.get(rule_name)
