
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
        "problem_smarts": "[NH2][c]",
        "candidates": [
            {"smiles": "C(=O)N", "name": "acetamide (acylated amine)",
             "rationale": "1차 방향족 아민을 아마이드로 아실화하여 N-hydroxylation 경로 자체를 차단"},
            {"smiles": "F", "name": "fluorine",
             "rationale": "반응성 아민을 제거하면서 전자끄는기로 고리 전자밀도 보정"},
        ],
    },
}

def get_replacement_candidates(rule_name: str) -> dict | None:
    """rule_name에 해당하는 치환 정보(SMARTS + 후보 리스트)를 반환. 없으면 None."""
    return REPLACEMENT_LIBRARY.get(rule_name)
