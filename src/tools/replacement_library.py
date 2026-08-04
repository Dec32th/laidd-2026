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
        "edit_method": "atom_edit",
        "problem_smarts": "[CX3H1](=O)",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "add_substituent", "param": "N",
             "target_idx_in_pattern": 0,
             "name": "amide",
             "rationale": "알데히드의 친전자성(단백질 부가물 형성 우려)을 제거하면서 "
                      "유사한 형태 유지. atom_edit 방식으로 재설계(기존 fragment-cut "
                      "은 회전 가능 결합으로 분리되지 않는 특수 맥락, 예: 폼아마이드형 "
                      "알데히드에서 조각화 실패)."},
            {"edit_type": "reduce_bond", "target_idx_pair_in_pattern": (0, 1),
             "name": "alcohol",
             "rationale": "가장 단순한 환원형 대체, 반응성 크게 감소. atom_edit 방식으로 "
                      "재설계(기존 fragment-cut 한계 해결)."},
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
        "problem_smarts": "[#6]S(=O)(=O)[OX2H1,OX1-]",
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
                      "아니라 약효 상실로 이어짐. 실제 도파민은 도파민 수용체 "
                      "D1(Ki 4.3-5.6 nM), D2(Ki 4.7-7.2 nM), D3(Ki 6.4-7.3 nM)에 "
                      "단자릿수 나노몰 수준의 강력한 작용제 친화도를 가짐(IUPHAR/BPS "
                      "Guide to PHARMACOLOGY 확인). || 인체의 COMT(catechol-O-"
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
    "thiol_1_dithiocarbamate": {
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
    "thiol_1_thiocarboxylate": {
        "edit_method": "atom_edit",
        "problem_smarts": "[SX1-]C(=O)",
        "target_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "replace_element", "param": 8, "name": "carboxylate (O replacing S)",
             "rationale": "티오카르복실산 음이온(R-C(=O)-S-)의 황을 산소로 대체하여 "
                          "카르복실산염(R-C(=O)-O-)으로 전환. 황 원자의 금속 킬레이팅 "
                          "및 친핵성 반응성을 제거함 (검증 필요)"},
        ],
    },
    "het-C-het_not_in_ring": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CX4]([OX2,SX2])([OX2,SX2])",
        "candidates": [
            {"edit_type": "remove_substituent",
             "center_idx_in_pattern": 0,
             "remove_idx_in_pattern": 1,
             "upgrade_bond_to_idx_in_pattern": 2,
             "name": "ketone/ester (one heteroatom substituent removed, C=O formed)",
             "rationale": "아세탈/케탈/오르토에스터(산소 2개) 또는 디티오아세탈(황 2개, "
                          "실제 철수약물 Probucol에서 확인) 등 탄소 하나에 헤테로원자 2개가 "
                          "붙은 구조는 가수분해/해리에 민감하여 반응성 카르보닐로 쉽게 "
                          "전환되며 대사 불안정성을 일으킴. 헤테로원자 하나를 제거하고 "
                          "남은 것을 카르보닐로 승격시켜, 가수분해로 어차피 도달할 안정한 "
                          "최종 형태로 미리 전환함 (검증 필요)"},
        ],
    },
    "cyclic_imide": {
        "edit_method": "atom_edit",
        "problem_smarts": "[C;R](=O)[N;R][C;R](=O)",
        "candidates": [
            {"edit_type": "cleave_bond", "cleave_pair_in_pattern": (2, 3),
             "name": "ring-opened amide (imide bond cleaved)",
             "rationale": "고리형 이미드(우레이드) 구조는 바르비투레이트류(페노바르비탈, "
                      "펜토바르비탈 등 다수 철수약물에서 실제 확인됨)와 탈리도마이드의 "
                      "잔여 글루타르이미드 고리에서 나타나며, 가수분해에 민감한 반응성 "
                      "구조임. 고리 내 아마이드 결합 하나를 끊어 개환함으로써 실제 "
                      "가수분해의 첫 단계를 근사함. 고리 구성원(R)만 매치하도록 제한하여, "
                      "개환 후 남은 사슬에 재적용되어 조각화되는 것을 방지함 "
                      "(ChEMBL 조회로 검증된 실제 철수약물 다수에서 발견, 검증 필요)"},
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
    "diketo_group": {
        "edit_method": "atom_edit",
        "problem_smarts": "C(=O)C(=O)",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "alpha-hydroxy ketone (reduced)",
             "rationale": "비시날 알파-디케톤(1,2-diketone)은 반응성이 높은 친전자체로 "
                          "단백질과 부가물을 형성할 수 있으며, 흡입 시 호흡기 독성을 "
                          "일으키는 것으로 알려진 디아세틸(버터향 첨가제) 사례가 대표적임. "
                          "카르보닐 하나를 환원하여 알파-하이드록시케톤(아실로인)으로 "
                          "전환, 케토-환원효소에 의한 실제 해독 경로와 유사한 방향으로 "
                          "반응성을 낮춤 (검증 필요)"},
        ],
    },
    "thioester": {
        "edit_method": "atom_edit",
        "problem_smarts": "[SX2](C(=O))",
        "target_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "replace_element", "param": 8, "name": "ester (O replacing S)",
             "rationale": "티오에스터의 황을 산소로 대체하여 일반 에스터로 전환. "
                          "티오에스터는 일반 에스터보다 가수분해 반응성이 높고 아실화 "
                          "능력이 강해 단백질 등과 부반응 우려가 있음 (검증 필요)"},
        ],
    },
    "N-nitroso": {
        "edit_method": "atom_edit",
        "problem_smarts": "[NX2;+0;!$(N(=O)[O-])]=[OX1;+0]",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "N-hydroxylamine (reduced)",
             "rationale": "N-니트로소 화합물(니트로사민)은 대사 활성화(알파-수산화)를 "
                          "거쳐 강력한 알킬화 발암물질을 생성하는 것으로 잘 알려짐 "
                          "(발사르탄, 라니티딘 등 실제 의약품 불순물 리콜 사례). "
                          "N=O를 환원하여 반응성을 낮춤 (검증 필요: 완전한 해독은 "
                          "탈니트로소화가 필요하며 이는 근사적 접근)"},
        ],
    },
    "hydrazine": {
        "edit_method": "atom_edit",
        "problem_smarts": "[NX3H2][NX3H1]",
        "center_idx_in_pattern": 1,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 0,
             "center_idx_in_pattern": 1,
             "name": "amide/amine (terminal N removed)",
             "rationale": "하이드라진/하이드라지드(R-NH-NH2)의 말단 질소를 제거하여 "
                          "단순 아민 또는 아마이드로 되돌림. 하이드라진류는 대사 시 "
                          "반응성 디아제늄 중간체를 형성해 유전독성을 일으킬 수 있는 "
                          "것으로 알려짐. 이는 azo_A(324) 환원 시 생성되는 하이드라진 "
                          "중간체의 잔여 위험을 추가로 낮추는 후속 규칙이기도 함 "
                          "(검증 필요)"},
        ],
    },
    "sulphate": {
        "edit_method": "atom_edit",
        "problem_smarts": "[OX2][SX4](=O)(=O)[OX1,OX2H]",
        "center_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 1,
             "center_idx_in_pattern": 0,
             "name": "alcohol (sulfate group removed)",
             "rationale": "알킬 설페이트 에스터(R-O-SO3-)는 대사되어 반응성 있는 "
                          "설페이트 이탈기를 통한 알킬화제로 작용할 수 있음(디메틸설페이트가 "
                          "강력한 발암/독성 물질로 잘 알려진 대표 사례). 설페이트기 전체를 "
                          "제거하여 원래의 알코올로 되돌림 (검증 필요)"},
        ],
    },
    "N_oxide": {
        "edit_method": "atom_edit",
        "problem_smarts": "[n+][O-]",
        "center_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 1,
             "center_idx_in_pattern": 0,
             "name": "pyridine (N-oxide removed)",
             "rationale": "방향족 N-옥사이드는 산화적 대사산물이자 반응성 중간체 "
                          "생성 경로의 일부일 수 있음. 산소를 제거하여 원래의 중성 "
                          "방향족 아민(피리딘 등)으로 환원, 자연 대사에서의 환원 "
                          "경로와 유사한 방향으로 반응성을 낮춤 (검증 필요)"},
        ],
    },
    "2-halo_pyridine": {
        "edit_method": "atom_edit",
        "problem_smarts": "n:c(-[Cl,Br,I])",
        "center_idx_in_pattern": 1,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 2,
             "center_idx_in_pattern": 1,
             "name": "pyridine (halogen removed)",
             "rationale": "피리딘 고리 질소에 인접한 위치의 할로겐(특히 불소/염소)은 "
                          "친핵성 방향족 치환(SNAr) 반응에 취약해, 체내 친핵체(글루타치온, "
                          "단백질 시스테인 등)와 반응할 수 있음. 할로겐을 제거하고 수소로 "
                          "대체하여 이 반응성 경로를 차단함 (검증 필요)"},
        ],
    },
    "disulphide": {
        "edit_method": "atom_edit",
        "problem_smarts": "[SX2][SX2]",
        "candidates": [
            {"edit_type": "cleave_bond", "cleave_pair_in_pattern": (0, 1),
             "name": "two thiols (bond cleaved)",
             "rationale": "[참고] 이황화결합(S-S)은 시스틴/단백질의 3차구조 형성에 "
                          "필수적인 정상 생체 구조이기도 하므로, 이 결합이 약물의 "
                          "구조 안정성이나 표적 결합에 관여하는 경우 본 치환이 "
                          "부적절할 수 있음. || 디티오카바메이트류(티우람 등) 농약/"
                          "살균제에서 흔한 반응성 이황화결합을 두 개의 티올로 분리, "
                          "산화·금속킬레이팅 반응성을 낮춤 (검증 필요)"},
        ],
    },
    "quinone_A(370)": {
        "edit_method": "atom_edit",
        "problem_smarts": "O=C1C=CC(=O)C=C1",
        "target_pairs_in_pattern": [(1, 0), (4, 5)],
        "ring_atoms_in_pattern": [1, 2, 3, 4, 6, 7],
        "ring_bonds_in_pattern": [(1, 2), (2, 3), (3, 4), (4, 6), (6, 7), (7, 1)],
        "candidates": [
            {"edit_type": "reduce_multi_bond", "name": "hydroquinone (reduced, re-aromatized)",
             "rationale": "파라벤조퀴논은 산화환원 사이클(redox cycling)을 통해 활성산소종(ROS)을 "
                          "생성하고 DNA/단백질과 직접 공유결합하는 대표적 반응성 구조. 체내 "
                          "NQO1(퀴논 환원효소) 효소가 실제로 수행하는 반응과 동일하게 두 카르보닐을 "
                          "환원하고 고리를 재방향족화하여 안정적인 하이드로퀴논으로 전환. 결과물이 "
                          "다시 hydroquinone 규칙에 해당할 수 있으며, 이 경우 반복 루프가 자동으로 "
                          "메톡시페놀 등 산화에 더 안정적인 형태로 한 단계 더 개선함 (검증 필요, "
                          "안트라퀴논 등 융합고리형은 미지원)"},
        ],
    },
    "quinone_A_anthraquinone": {
        "edit_method": "atom_edit",
        "problem_smarts": "O=C1c2ccccc2C(=O)c2ccccc21",
        "target_pairs_in_pattern": [(1, 0), (8, 9)],
        "ring_atoms_in_pattern": [1, 2, 3, 4, 5, 6, 7, 8],
        "ring_bonds_in_pattern": [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 1)],
        "candidates": [
            {"edit_type": "reduce_multi_bond", "name": "anthrahydroquinone (reduced, re-aromatized)",
             "rationale": "안트라퀴논은 벤조퀴논과 동일한 산화환원 사이클링(redox cycling) 메커니즘을 "
                          "가지되, 두 벤젠 고리에 의해 안정화되어 항암제(독소루비신 등) 및 염료에서도 "
                          "흔히 쓰이는 골격임. 두 카르보닐을 동시에 환원하고 중앙 고리를 재방향족화하여 "
                          "안트라하이드로퀴논으로 전환, 산화환원 사이클링 능력을 제거함. 결과물이 "
                          "hydroquinone 규칙에 해당할 수 있어 반복 루프가 자동으로 추가 개선 가능 "
                          "(Murcko scaffold 분석으로 발견, 검증 필요)"},
        ],
    },
    "quinone_diimine": {
        "edit_method": "atom_edit",
        "problem_smarts": "N=C1C=CC(=N)C=C1",
        "target_pairs_in_pattern": [(1, 0), (4, 5)],
        "ring_atoms_in_pattern": [1, 2, 3, 4, 6, 7],
        "ring_bonds_in_pattern": [(1, 2), (2, 3), (3, 4), (4, 6), (6, 7), (7, 1)],
        "candidates": [
            {"edit_type": "reduce_multi_bond", "name": "phenylenediamine (reduced, re-aromatized)",
             "rationale": "퀴논디이민(quinone diimine)은 벤조퀴논의 산소가 이민으로 치환된 유사체로, "
                          "동일한 산화환원 사이클링 메커니즘을 가지며 헤어염료 성분(파라페닐렌디아민 "
                          "산화형) 등에서 피부 알레르기 및 접촉성 피부염을 유발하는 것으로 알려짐. 두 "
                          "이민을 동시에 환원하고 고리를 재방향족화하여 페닐렌디아민(원래의 안정한 "
                          "환원형)으로 전환 (Murcko scaffold 분석으로 발견, 검증 필요)"},
        ],
    },
    "isocyanate": {
        "edit_method": "atom_edit",
        "problem_smarts": "[NX2]=[CX2]=[OX1]",
        "center_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 1,
             "center_idx_in_pattern": 0,
             "name": "amine (NCO hydrolyzed)",
             "rationale": "이소시아네이트(R-N=C=O)는 매우 반응성이 높은 친전자체로, "
                          "단백질/아미노기와 쉽게 부가반응을 일으켜 직업성 천식·과민증을 "
                          "유발하는 것으로 잘 알려짐(TDI, MDI 등 산업용 이소시아네이트 "
                          "사례). 체내/환경에서 실제로 일어나는 가수분해 경로(R-NCO + H2O "
                          "-> R-NH2 + CO2)와 동일하게 카르보닐 탄소와 산소를 제거하고 "
                          "질소만 남겨 아민으로 전환 (검증 필요)"},
        ],
    },
    "triple_bond": {
        "problem_smarts": "C#C",
        "edit_method": "atom_edit",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "alkene (partially reduced)",
             "rationale": "말단 알카인(삼중결합)은 CYP450 효소에 의해 기계기반 억제"
                          "(mechanism-based inhibition) 경로로 대사되며, 반응성 케텐/"
                          "에폭사이드 중간체를 형성해 효소를 비가역적으로 불활성화할 "
                          "수 있음(에티닐에스트라디올 등에서 알려진 메커니즘). 삼중결합을 "
                          "이중결합으로 환원하여 반응성을 낮춤 (검증 필요, 완전 포화가 "
                          "아닌 부분 환원)"},
        ],
    },
    "stilbene": {
        "problem_smarts": "c-[CX3]=[CX3]-c",
        "edit_method": "atom_edit",
        "target_idx_pair_in_pattern": (1, 2),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "diarylethane (reduced)",
             "rationale": "스틸벤 구조(두 방향족 고리를 잇는 C=C)는 디에틸스틸베스트롤"
                          "(DES)처럼 내분비교란 및 대사 산화를 통한 반응성 중간체 형성이 "
                          "알려진 골격. 이중결합을 환원하여 평면성을 낮추고 대사 반응성을 "
                          "완화함 (검증 필요, 에스트로겐 수용체 결합에 필요한 형태 자체를 "
                          "훼손할 수 있어 신중한 해석 필요)"},
        ],
    },
    "beta-keto/anhydride": {
        "edit_method": "atom_edit",
        "problem_smarts": "C(=O)OC(=O)",
        "center_idx_in_pattern": 2,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 3,
             "center_idx_in_pattern": 2,
             "name": "carboxylic acid (anhydride hydrolyzed)",
             "rationale": "산 무수물(R-C(=O)-O-C(=O)-R')은 강한 아실화제로 단백질 아미노산 "
                          "잔기와 쉽게 반응하며, 수용액 환경에서 자발적으로 가수분해되어 "
                          "두 개의 카르복실산으로 분해되는 것이 자연스러운 무독화 경로임. "
                          "한쪽 아실기를 제거하여 이 가수분해 최종형(카르복실산)으로 직접 "
                          "전환 (검증 필요). ※ 대안 후보(무수물->아마이드/이미드 bioisostere) "
                          "는 문헌 확인 후 추가 예정"},
        ],
    },
    "phthalimide": {
        "edit_method": "atom_edit",
        "problem_smarts": "O=C1c2ccccc2C(=O)N1[#6]",
        "center_idx_in_pattern": 10,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 1,
             "center_idx_in_pattern": 10,
             "name": "primary amine (imide hydrolyzed)",
             "rationale": "프탈이미드(고리형 이미드)는 탈리도마이드 등에서 알려진 골격으로, "
                          "체내에서 가수분해되어 원래의 1차 아민과 프탈산으로 분해되는 것이 "
                          "자연스러운 대사 경로임. 이 가수분해 용이성 자체가 대사 불안정성/"
                          "반응성 우려의 근거이며, 고리 전체를 제거하여 이 가수분해 최종형인 "
                          "1차 아민으로 직접 전환 (검증 필요)"},
        ],
    },
    "hydroxamic_acid": {
        "edit_method": "atom_edit",
        "problem_smarts": "C(=O)N[OX2H1]",
        "center_idx_in_pattern": 2,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 3,
             "center_idx_in_pattern": 2,
             "name": "amide (N-hydroxyl removed)",
             "rationale": "[참고] 하이드록삼산(R-C(=O)-NH-OH)은 보리노스타트, 파노비노스타트 "
                          "등 HDAC 억제제에서 아연 킬레이션을 통한 핵심 약효 작용기로 쓰이므로, "
                          "이 계열에는 본 치환이 약효 상실로 이어질 수 있음. || 하이드록삼산은 "
                          "로센 재배열(Lossen rearrangement)을 통해 반응성 이소시아네이트로 "
                          "전환될 수 있는 잠재적 위험이 있음. N-하이드록실기를 제거해 단순 "
                          "아마이드로 전환, 이 재배열 경로를 차단함 (검증 필요)"},
        ],
    },
    "Aliphatic_long_chain": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CH2][CH2][CH2][CH2]",
        "candidates": [
            {"edit_type": "insert_atom", "insert_pair_in_pattern": (1, 2), "param": 8,
             "name": "ether-inserted chain (O in middle)",
             "rationale": "..."},  # 기존 그대로 유지
            {"edit_type": "insert_atom_multi_chain",
             "chain_start_idx_in_pattern": 0,
             "name": "multi-ether chain (multiple O inserted for long chains)",
             "rationale": "매우 긴 지방족 사슬(수 회 반복이 필요한 경우)에 대해, 4탄소 "
                          "간격마다 산소를 동시에 여러 개 삽입하여 한 번에 극성을 분산시킴 "
                          "(검증 필요, 긴 사슬 전용)"},
        ],
    },
    "isolated_alkene": {
        "edit_method": "atom_edit",
        "problem_smarts": "[CX3H1,CX3H0;!$([CX3]=[CX3]c)]=[CX3;!$([CX3]=[CX3]c)]",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "saturated (C-C single bond)",
             "rationale": "고립된 지방족 알켄(방향족·카르보닐과 공액되지 않은 단순 C=C)은 "
                          "산화적 대사(에폭시드 형성 등)를 거쳐 반응성 중간체를 생성할 "
                          "가능성이 있는 구조 경고임. 이중결합을 단일결합으로 환원해 이 "
                          "산화 경로를 차단함 (검증 필요)"},
        ],
    },
    "quaternary_nitrogen_1": {
        "edit_method": "atom_edit",
        "problem_smarts": "[#6][n+]1ccccc1",
        "center_idx_in_pattern": 1,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 0,
             "center_idx_in_pattern": 1,
             "allow_counterion": True,
             "allow_aromatic_zero_h": True,
             "name": "pyridine (N-alkyl removed)",
             "rationale": "N-알킬피리디늄(방향족 4차 질소)은 영구적 양전하를 띠어 세포막 "
                          "투과성이 떨어지고, 파라쿼트 등 일부 사례에서 미토콘드리아 "
                          "독성/신경독성과 연관됨. N-알킬 사슬을 제거해 중성 피리딘으로 "
                          "복원함 (검증 필요)"},
        ],
    },
    "quaternary_nitrogen_2": {
        "edit_method": "atom_edit",
        "problem_smarts": "[#6][CH2][N+]([#6])([#6])[#6]",
        "center_idx_in_pattern": 2,
        "candidates": [
            {"edit_type": "remove_atom",
             "remove_idx_in_pattern": 1,
             "center_idx_in_pattern": 2,
             "name": "tertiary amine (one alkyl removed)",
             "rationale": "비방향족 4차 암모늄(영구적 양전하)은 신경근 차단제(예: "
                          "석시닐콜린류)에서 보이는 것처럼 막 투과성 저하 및 특정 이온"
                          "채널/수용체와의 비특이적 상호작용 우려가 있음. 알킬기 하나를 "
                          "제거해 중성 3차 아민으로 복원함 (검증 필요)"},
        ],
    },
    "phenol_ester": {
        "edit_method": "atom_edit",
        "problem_smarts": "c[OX2]C(=O)",
        "candidates": [
            {"edit_type": "cleave_bond", "cleave_pair_in_pattern": (0, 1),
             "name": "phenol + carboxylic acid (ester cleaved)",
             "rationale": "페놀 에스터(아릴-O-C(=O)-)는 일반 지방족 에스터보다 가수분해에 "
                          "민감하고, 방출되는 페놀이 추가로 반응성 퀴논으로 산화될 수 있는 "
                          "이중 우려가 있는 구조임. 에스터 결합을 끊어 페놀과 카르복실산으로 "
                          "분리, 가수분해로 어차피 도달할 안정한 최종 형태로 전환 (검증 필요)"},
        ],
    },
    "phosphor": {
        "edit_method": "atom_edit",
        "problem_smarts": "[OX2][PX4](=[OX1])([OX2])[OX2]",
        "candidates": [
            {"edit_type": "cleave_bond", "cleave_pair_in_pattern": (1, 4),
             "name": "diester + phenol/alcohol (one ester bond cleaved)",
             "rationale": "유기인산 트리에스터(트리아릴/트리알킬 포스페이트)는 아세틸콜린"
                          "에스터라제(AChE) 억제를 통한 신경독성 메커니즘이 잘 알려진 "
                          "구조로(유기인계 살충제·신경작용제의 공통 골격), 다중 에스터 "
                          "결합이 반응성/생체이용률에 기여함. 에스터 결합 하나를 가수분해로 "
                          "끊어 반응성을 낮춤 (검증 필요, 인 원자에 남은 나머지 에스터는 "
                          "추가 규칙 필요 가능)"},
        ],
    },
}

def get_replacement_candidates(rule_name: str) -> dict | None:
    """rule_name에 해당하는 치환 정보(SMARTS + 후보 리스트)를 반환. 없으면 None."""
    return REPLACEMENT_LIBRARY.get(rule_name)
