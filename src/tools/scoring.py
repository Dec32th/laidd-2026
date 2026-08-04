"""치환 후보에 대한 다목적 종합 스코어.
Toxicity/SA Score/QED/Lipinski/PAINS는 즉시 계산, Docking은 선례
라이브러리에 이미 검증된 값이 있을 때만 반영(실시간 도킹은 시간상
불가능 - 오늘 세션에서 1건당 수 분~수십 분 소요, 배치 처리에 부적합).
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, FilterCatalog


def _lipinski_violations(mol):
    violations = 0
    if Descriptors.MolWt(mol) > 500: violations += 1
    if Descriptors.MolLogP(mol) > 5: violations += 1
    if Descriptors.NumHDonors(mol) > 5: violations += 1
    if Descriptors.NumHAcceptors(mol) > 10: violations += 1
    return violations


def _pains_pass(mol, catalog):
    return not catalog.HasMatch(mol)


def compute_multi_objective_score(original_smiles, fixed_smiles, rule_name,
                                    tox_delta=None, sascorer_module=None,
                                    precedent_docking_delta=None,
                                    weights=None):
    """0~1 범위로 정규화한 항목별 점수와 가중합을 반환.
    tox_delta: 외부에서 baseline 모델로 계산한 독성 예측값 변화(음수=개선),
               없으면 None으로 두고 해당 항목 제외.
    sascorer_module: sascorer 모듈(외부에서 import해서 전달, 순환import 방지).
    precedent_docking_delta: 선례 라이브러리에 해당 규칙의 도킹 kcal/mol
               변화값이 있으면 전달(음수=결합강화), 없으면 None.
    weights: 항목별 가중치 딕셔너리, 기본값 아래 참고."""
    default_weights = {"toxicity": 0.30, "docking": 0.20, "sa": 0.15,
                        "qed": 0.15, "lipinski": 0.10, "pains": 0.10}
    w = weights or default_weights

    mol_o = Chem.MolFromSmiles(original_smiles)
    mol_f = Chem.MolFromSmiles(fixed_smiles)
    if mol_o is None or mol_f is None:
        return None

    scores = {}
    used_weight = 0.0

    if tox_delta is not None:
        scores["toxicity"] = max(0.0, min(1.0, 0.5 - tox_delta))
        used_weight += w["toxicity"]

    if precedent_docking_delta is not None:
        scores["docking"] = max(0.0, min(1.0, 0.5 - precedent_docking_delta / 2))
        used_weight += w["docking"]

    if sascorer_module is not None:
        sa_o = sascorer_module.calculateScore(mol_o)
        sa_f = sascorer_module.calculateScore(mol_f)
        scores["sa"] = max(0.0, min(1.0, 1 - (sa_f - sa_o) / 5))
        used_weight += w["sa"]

    qed_o, qed_f = QED.qed(mol_o), QED.qed(mol_f)
    scores["qed"] = max(0.0, min(1.0, 0.5 + (qed_f - qed_o)))
    used_weight += w["qed"]

    lip_o = _lipinski_violations(mol_o)
    lip_f = _lipinski_violations(mol_f)
    scores["lipinski"] = max(0.0, min(1.0, 0.5 + (lip_o - lip_f) * 0.25))
    used_weight += w["lipinski"]

    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    scores["pains"] = 1.0 if _pains_pass(mol_f, catalog) else 0.0
    used_weight += w["pains"]

    weighted_sum = sum(scores[k] * w[k] for k in scores)
    composite = weighted_sum / used_weight if used_weight > 0 else None

    return {
        "component_scores": scores,
        "weights_used": {k: w[k] for k in scores},
        "composite_score": composite,
        "note": "docking/toxicity는 선례·외부계산 존재 시에만 반영, 부재 시 나머지 항목으로 정규화",
    }
