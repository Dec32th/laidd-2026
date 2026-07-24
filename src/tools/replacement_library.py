
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
    "michael_acceptor": {
        "problem_smarts": "C=CC(=O)",
        "candidates": [
            {"smiles": "CCC(=O)", "name": "saturated ketone",
             "rationale": "이중결합을 제거해 단백질 친전자성 부가반응(covalent binding) 위험 제거"},
        ],
    },
    "thiourea": {
        "problem_smarts": "NC(=S)N",
        "candidates": [
            {"smiles": "NC(=O)N", "name": "urea",
             "rationale": "황 원자를 산소로 대체, 유사한 형태를 유지하면서 반응성/대사 우려 감소"},
        ],
    },
    "acyl_halide": {
        "problem_smarts": "C(=O)[F,Cl,Br,I]",
        "candidates": [
            {"smiles": "C(=O)N", "name": "amide",
             "rationale": "고반응성 아실할라이드를 안정적인 아마이드로 대체"},
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
    "phenol": {
        "problem_smarts": "[OX2H]",
        "candidates": [
            {"smiles": "Cl", "name": "chlorine",
             "rationale": "MMPA 데이터 기반 발견: Tox21에서 페놀->염소 치환 시 NR-ER, "
                          "NR-ER-LBD, SR-ARE 3개 assay 동시 개선 관찰됨. 페놀의 산화적 대사 "
                          "(quinone 형성 등) 경로를 차단하는 것으로 추정",
             "source": "data-driven (Tox21 MMP analysis, charge-matched, Δmw<50)"},
            {"smiles": "C", "name": "methyl",
             "rationale": "MMPA 데이터 기반 발견: Tox21에서 페놀->메틸 치환 시 NR-AR "
                          "assay 개선 관찰됨. 히드록실기 제거로 산화 취약성 감소 추정",
             "source": "data-driven (Tox21 MMP analysis, charge-matched, Δmw<50)"},
        ],
    },
    "amide": {
        "problem_smarts": "[NX3H1][CX3](=O)[#6]",
        "candidates": [
            {"smiles": "NC(=O)N", "name": "urea",
             "rationale": "MMPA 데이터 기반 발견: Tox21에서 N-메틸아마이드->우레아 치환 시 "
                          "SR-ARE assay 개선 관찰됨. 기존 thiourea->urea 치환과 같은 "
                          "계열(우레아 활용)로 일관성 있음",
             "source": "data-driven (Tox21 MMP analysis, charge-matched, Δmw<50)"},
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
        ],
    },
    "imine_1": {
        "problem_smarts": "C=N[OX2H1]",
        "candidates": [
            {"smiles": "CN", "name": "amine (reduced)",
             "rationale": "옥심의 C=N 결합을 환원하여, 가수분해 시 원래의 반응성 "
                          "카르보닐(알데히드/케톤)로 되돌아갈 수 있는 대사 불안정 "
                          "경로를 제거함"},
        ],
    },
}

def get_replacement_candidates(rule_name: str) -> dict | None:
    """rule_name에 해당하는 치환 정보(SMARTS + 후보 리스트)를 반환. 없으면 None."""
    return REPLACEMENT_LIBRARY.get(rule_name)
