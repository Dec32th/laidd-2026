"""치환 전후 활성 보존 가능성을 정성/정량으로 평가하는 지표 모음.
2D 연결성(Tanimoto), 3D 형태(회전반경), QED/LogP/합성용이성 변화를
종합해 판정한다. sascorer는 외부 저작물이라 파라미터로 주입받는다
(scoring.py와 동일한 패턴)."""

from rdkit import Chem
from rdkit.Chem import Descriptors, Descriptors3D, DataStructs, QED, AllChem
from models.tox_baseline import smiles_to_fp_bitvect


def compute_activity_preservation_metrics(original_smiles, fixed_smiles, sascorer_module):
    mol_o = Chem.MolFromSmiles(original_smiles)
    mol_f = Chem.MolFromSmiles(fixed_smiles)
    if mol_o is None or mol_f is None:
        return None

    fp_o = smiles_to_fp_bitvect(original_smiles)
    fp_f = smiles_to_fp_bitvect(fixed_smiles)
    tanimoto = DataStructs.TanimotoSimilarity(fp_o, fp_f)

    qed_o, qed_f = QED.qed(mol_o), QED.qed(mol_f)
    logp_o, logp_f = Descriptors.MolLogP(mol_o), Descriptors.MolLogP(mol_f)
    sa_o = sascorer_module.calculateScore(mol_o)
    sa_f = sascorer_module.calculateScore(mol_f)

    def get_3d(mol):
        m = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(m)
        return m

    m3d_o, m3d_f = get_3d(mol_o), get_3d(mol_f)
    shape_available = m3d_o is not None and m3d_f is not None

    result = {"tanimoto": tanimoto, "delta_qed": qed_f - qed_o, "delta_logp": logp_f - logp_o,
              "delta_sa_score": sa_f - sa_o, "shape_available": shape_available}
    if shape_available:
        rog_o = Descriptors3D.RadiusOfGyration(m3d_o)
        rog_f = Descriptors3D.RadiusOfGyration(m3d_f)
        result["delta_rog_pct"] = (rog_f - rog_o) / rog_o * 100 if rog_o != 0 else None
    return result


def classify_activity_risk_v3(metrics):
    if metrics is None:
        return {"verdict": "판정 불가", "details": []}
    details, warnings = [], []

    conn_ok = metrics["tanimoto"] >= 0.5
    details.append(f"2D 연결성: {'유사' if conn_ok else '상이'} (Tanimoto {metrics['tanimoto']:.3f})")

    shape_ok = None
    if metrics["shape_available"] and metrics.get("delta_rog_pct") is not None:
        shape_ok = abs(metrics["delta_rog_pct"]) < 15
        details.append(f"3D 형태: {'보존' if shape_ok else '변화'} (회전반경 {metrics['delta_rog_pct']:+.1f}%)")
    else:
        details.append("3D 형태: 계산 불가")

    qed_ok = abs(metrics["delta_qed"]) < 0.1
    details.append(f"약물유사성(QED): {'유지' if qed_ok else '변화'} ({metrics['delta_qed']:+.3f})")
    if not qed_ok:
        warnings.append("QED 변화")

    logp_ok = abs(metrics["delta_logp"]) < 1.0
    details.append(f"소수성(LogP): {'유지' if logp_ok else '변화'} ({metrics['delta_logp']:+.3f})")
    if not logp_ok:
        warnings.append("LogP 변화")

    sa_ok = metrics["delta_sa_score"] < 0.5
    details.append(f"합성용이성(SA): {'유지/개선' if sa_ok else '악화'} ({metrics['delta_sa_score']:+.3f})")
    if not sa_ok:
        warnings.append("합성난이도 증가")

    if shape_ok is None:
        verdict = "3D 형태 계산 불가 — 2D 지표만으로 판단, 신뢰도 낮음"
    elif shape_ok:
        verdict = ("구조·형태 모두 보존 — 활성 유지 가능성 높음" if conn_ok else
                   "2D 연결성은 크게 바뀌었으나 3D 형태는 보존됨 (bioisostere 가능성) — 활성 유지 기대")
    else:
        verdict = "3D 형태 자체가 크게 변화 — 표적 결합 형태 훼손 우려, 사람 검토 필요"

    if warnings:
        verdict += f" [보조 경고: {', '.join(warnings)}]"

    return {"verdict": verdict, "details": details, "shape_ok": shape_ok, "warnings": warnings}
