import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import train_test_split

TOX21_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"

_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _is_valid_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None

def _smiles_to_ecfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = _generator.GetFingerprintAsNumPy(mol)
    return np.array(fp)

def load_tox21_clean(test_size=0.3, valid_ratio=0.5, random_state=42):
    """Tox21 데이터셋을 로드하고 train/valid/test로 분할.

    SMILES 파싱 실패 분자는 제외하고, ECFP 지문(fingerprint)과
    12개 assay 레이블(y), 결측 마스크(w)를 함께 준비한다.
    test_size는 전체 중 (valid+test)에 할당할 비율, valid_ratio는
    그 중 test로 다시 나눌 비율(기본: valid/test 각 절반씩).

    Returns:
        dict: smiles_train/valid/test, X_*, y_*, w_*, invalid_smiles
        등을 담은 딕셔너리.
    """
    df = pd.read_csv(TOX21_URL)
    df['valid'] = df['smiles'].apply(_is_valid_smiles)

    n_total, n_valid = len(df), df['valid'].sum()
    print(f"전체: {n_total}개, 파싱 성공: {n_valid}개, 파싱 실패(제외): {n_total - n_valid}개")

    invalid_smiles = df[~df['valid']]['smiles'].tolist()
    df_clean = df[df['valid']].reset_index(drop=True)

    task_cols = [c for c in df.columns if c not in ['smiles', 'mol_id', 'valid']]

    y = df_clean[task_cols].fillna(0).values.astype(np.float32)
    w = (~df_clean[task_cols].isna()).values.astype(np.float32)
    X = np.stack(df_clean['smiles'].apply(_smiles_to_ecfp).values)
    smiles_arr = df_clean['smiles'].values

    indices = np.arange(len(X))
    train_idx, temp_idx = train_test_split(indices, test_size=test_size, random_state=random_state)
    valid_idx, test_idx = train_test_split(temp_idx, test_size=valid_ratio, random_state=random_state)

    return {
        'X_train': X[train_idx], 'y_train': y[train_idx], 'w_train': w[train_idx],
        'smiles_train': smiles_arr[train_idx],
        'X_valid': X[valid_idx], 'y_valid': y[valid_idx], 'w_valid': w[valid_idx],
        'smiles_valid': smiles_arr[valid_idx],
        'X_test': X[test_idx], 'y_test': y[test_idx], 'w_test': w[test_idx],
        'smiles_test': smiles_arr[test_idx],
        'task_cols': task_cols,
        'invalid_smiles': invalid_smiles,
    }
