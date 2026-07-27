
REPLACEMENT_LIBRARY = {
    "nitro_group": {
        "problem_smarts": "[N+](=O)[O-]",
        "candidates": [
            {"smiles": "N", "name": "primary amine",
             "rationale": "극성을 유지하면서 니트로기의 환원성 대사 중간체 생성 경로를 제거함"},
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
             "rationale": "알파,베타-불포화 카르보닐의 C=C 이중결합을 환원하여 "
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
             "rationale": "이탈기를 제거해 알킬화 반응성을 없앰, 극성은 유사하게 유지"},
            {"smiles": "F", "name": "fluorine",
             "rationale": "할로겐을 유지하되 C-F 결합은 강해 이탈기로 작용하지 않음, 입체적 크기도 유사"},
        ],
    },
    "aniline": {
        "problem_smarts": "[NH2]",
        "candidates": [
            {"smiles": "C(=O)N", "name": "acetamide (acylated amine)",
             "rationale": "1차 방향족 아민을 아마이드로 아실화하여 N-hydroxylation 경로 자체를 차단"},
            {"smiles": "F", "name": "fluorine",
             "rationale": "반응성 아민을 제거하면서 전자끄는기로 고리 전자밀도 보정"},
        ],
    },
    "Sulfonic_acid_2": {
        "problem_smarts": "S(=O)(=O)[OX2H1,OX1-]",
        "candidates": [
            {"smiles": "S(=O)(=O)N", "name": "sulfonamide",
             "rationale": "생리적 pH에서 이온화 정도(전하)를 크게 낮춰 세포막 투과성을 "
                          "개선함. 설폰산은 대부분 음이온 상태로 존재해 경구 흡수가 "
                          "저해되는 경우가 많으나, 설폰아마이드는 유사한 골격을 유지하면서도 "
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
        "problem_smarts": "C=N",
        "target_idx_pair_in_pattern": (0, 1),
        "candidates": [
            {"edit_type": "reduce_bond", "name": "amine (reduced)",
             "rationale": "일반 이민(C=N-R)을 환원하여 가수분해 시 반응성 카르보닐로 "
                          "되돌아갈 수 있는 대사 불안정 경로를 제거함. 옥심 특유의 "
                          "메커니즘보다는 근거가 다소 약하며, 하위 구조별 개별 검증 필요"},
        ],
    },
    "catechol": {
        "edit_method": "atom_edit",
        "problem_smarts": "[OX2H;$(Oc1ccccc1O)]",
        "target_idx_in_pattern": 0,
        "candidates": [
            {"edit_type": "add_substituent", "param": "C", "name": "methoxy",
             "rationale": "인체의 COMT(catechol-O-methyltransferase) 효소가 카테콜을 "
                          "메톡시페놀로 메틸화하여 해독하는 생리적 경로와 동일한 원리. "
                          "오르토-퀴논으로의 산화 경로를 차단하여 세포독성/유전독성 우려를 "
                          "낮춤 (학생 확인 예정: ScienceDirect catechol overview, "
                          "PMC6643002 등 참고)"},
        ],
    },
    "Thiocarbonyl_group": {
        "edit_method": "atom_edit",
        "problem_smarts": "[#6]=[#16]",
        "target_idx_in_pattern": 1,
        "candidates": [
            {"edit_type": "replace_element", "param": 8, "name": "carbonyl (O replacing S)",
             "rationale": "황을 산소로 대체(티오카르보닐->카르보닐)하는 것은 흔한 "
                          "bioisostere 전략으로, 갑상선 기능 저해 등 황 함유 작용기 "
                          "특유의 대사/독성 우려를 낮춤 (검증 필요, thiourea->urea "
                          "치환 논리와 동일 계열)"},
        ],
    },
    "aniline_ring_bcp": {
        "edit_method": "atom_edit",
        "problem_smarts": "[NH2]c1ccc([#6,#7,#8,#16])cc1",
        "ring_atom_indices_in_pattern": [1, 2, 3, 4, 6, 7],
        "anchor_indices_in_pattern": (0, 5),
        "candidates": [
            {"edit_type": "replace_ring", "param": "[*:1]C12CC(C1)(C2)[*:2]",
             "name": "BCP (bicyclo[1.1.1]pentane)",
             "rationale": "para-이치환 아닐린의 방향족 벤젠 고리를 포화 bicyclic "
                          "탄소골격(BCP)으로 교체함. 방향족성 제거로 aniline reactive "
                          "metabolite(RM) 형성 및 CYP-inhibition을 감소시켜, 퀴논이민 "
                          "생성 경로를 차단하고 특이체질 약물 부작용(IADR) 위험을 낮춤 "
                          "(문헌 근거, 학생 제공). 벤젠과의 공간적 유사성, Fsp3 증가, "
                          "실제 성공 사례가 많아 우선 채택함(BCO/NB/CUB는 근거 부족으로 보류)"},
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
}

def get_replacement_candidates(rule_name: str) -> dict | None:
    """rule_name에 해당하는 치환 정보(SMARTS + 후보 리스트)를 반환. 없으면 None."""
    return REPLACEMENT_LIBRARY.get(rule_name)
