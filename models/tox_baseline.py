"""Tox21 baseline 독성 예측 모델. 매 세션 직접 학습해서 사용(random_state
고정으로 재현성 확보). debate 로직의 tox_delta 콜백에 연결하기 위한 래퍼."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier

_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def smiles_to_ecfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return _generator.GetFingerprintAsNumPy(mol) if mol else None


def train_tox21_baseline(data):
    """data: load_tox21_clean() 반환 딕셔너리.
    반환: {"classifiers": dict, "task_cols": list} — 이 딕셔너리 자체가
    학습된 모델 전체이며, 세션 내내 이 변수 하나만 들고 다니면 됨."""
    X_train, y_train, w_train = data['X_train'], data['y_train'], data['w_train']
    task_cols = data['task_cols']
    classifiers = {}
    for i, task in enumerate(task_cols):
        train_mask = w_train[:, i] == 1
        clf = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)
        clf.fit(X_train[train_mask], y_train[train_mask, i])
        classifiers[task] = clf
    return {"classifiers": classifiers, "task_cols": task_cols}


def predict_tox21_avg(smiles, model):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = smiles_to_ecfp(smiles).reshape(1, -1)
    classifiers, task_cols = model["classifiers"], model["task_cols"]
    return float(np.mean([classifiers[t].predict_proba(fp)[0][1] for t in task_cols]))


def make_tox_predictor(model):
    """set_tox_predictor()에 바로 넘길 수 있는 콜백 생성.
    반환값은 (fixed 예측 - original 예측): 음수면 독성이 줄어든 것."""
    def _predict(original_smiles, fixed_smiles, rule_name):
        p_o = predict_tox21_avg(original_smiles, model)
        p_f = predict_tox21_avg(fixed_smiles, model)
        if p_o is None or p_f is None:
            return None
        return p_f - p_o
    return _predict
