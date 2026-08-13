"""도킹 자동화 파이프라인.
표적별 HETATM 처리(보조인자 유지 여부)는 반드시 사람이 먼저 확인해야
하므로, 새 표적을 추가할 때는 DOCKING_TARGETS에 keep_hetatm_codes를
명시적으로 등록하는 과정을 거친다(완전 자동화 금지 지점)."""

import os
import json
import re
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

VINA_BIN = os.path.join(os.getcwd(), "vina_bin")

DOCKING_TARGETS = {
    "catechol": {
        "target_name": "COMT", "pdb_id": "1VID",
        "ligand_code": "DNC", "keep_hetatm_codes": ["MG"],
        "box_size": [20, 20, 20],
    },
    "Michael_acceptor_1": {
        "target_name": "EGFR", "pdb_id": "6JX4",
        "ligand_code": "YY3", "keep_hetatm_codes": [],
        "box_size": [20, 20, 20],
    },
    "hydroquinone": {
        "target_name": "NQO1", "pdb_id": None,
        "ligand_code": None, "keep_hetatm_codes": ["FAD"],
        "box_size": [20, 20, 20],
    },
    "quinone_A(370)": {
        "target_name": "NQO1", "pdb_id": None,
        "ligand_code": None, "keep_hetatm_codes": ["FAD"],
        "box_size": [20, 20, 20],
    },
}

_CACHE_PATH = "outputs/docking_cache.json"


def _load_cache():
    if os.path.exists(_CACHE_PATH):
        with open(_CACHE_PATH) as f:
            return json.load(f)
    return {}


def _save_cache(cache):
    os.makedirs("outputs", exist_ok=True)
    with open(_CACHE_PATH, "w") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def prepare_ligand_pdbqt(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, f"{filename}.pdb")
    subprocess.run(["obabel", f"{filename}.pdb", "-O", f"{filename}.pdbqt"], capture_output=True)
    out_path = f"{filename}.pdbqt"
    return out_path if os.path.exists(out_path) else None


def run_docking_cli(ligand_pdbqt, receptor_pdbqt, out_prefix, box_center, box_size, exhaustiveness=4):
    log_path = f"{out_prefix}_log.txt"
    cmd = [
        VINA_BIN, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
        "--center_x", str(box_center[0]), "--center_y", str(box_center[1]), "--center_z", str(box_center[2]),
        "--size_x", str(box_size[0]), "--size_y", str(box_size[1]), "--size_z", str(box_size[2]),
        "--exhaustiveness", str(exhaustiveness), "--out", f"{out_prefix}_out.pdbqt",
    ]
    try:
        with open(log_path, "w") as log_f:
            subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, timeout=180)
    except subprocess.TimeoutExpired:
        return None

    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        log = f.read()
    match = re.search(r"^\s*1\s+(-?\d+\.\d+)", log, re.MULTILINE)
    return float(match.group(1)) if match else None


def prepare_receptor(pdb_id, keep_hetatm_codes):
    """수용체 준비. keep_hetatm_codes는 반드시 사람이 그 표적의 HETATM
    목록(prepare_receptor 호출 전 raw PDB를 먼저 열어 확인)을 보고
    직접 지정해야 하며, 빈 리스트라도 명시적으로 넘겨야 한다."""
    os.makedirs("targets", exist_ok=True)
    raw_path = f"targets/{pdb_id}.pdb"
    if not os.path.exists(raw_path):
        subprocess.run(["wget", "-q", f"https://files.rcsb.org/download/{pdb_id}.pdb", "-O", raw_path])

    with open(raw_path) as f:
        lines = f.readlines()

    keep_set = set(keep_hetatm_codes)
    clean_lines = [l for l in lines if l.startswith(("ATOM", "TER", "END"))
                   or (l.startswith("HETATM") and l[17:20].strip() in keep_set)]
    clean_path = f"targets/{pdb_id}_clean.pdb"
    with open(clean_path, "w") as f:
        f.writelines(clean_lines)

    receptor_pdbqt = f"targets/{pdb_id}_receptor.pdbqt"
    if not os.path.exists(receptor_pdbqt):
        subprocess.run(["obabel", clean_path, "-O", receptor_pdbqt, "-xr"], capture_output=True)
    return receptor_pdbqt, lines


def inspect_hetatm(pdb_id):
    """새 표적을 DOCKING_TARGETS에 등록하기 전, HETATM 목록을 먼저
    확인하는 용도. 자동 결정 없이 사람이 보고 keep_hetatm_codes를
    정하도록 정보만 제공한다."""
    os.makedirs("targets", exist_ok=True)
    raw_path = f"targets/{pdb_id}.pdb"
    if not os.path.exists(raw_path):
        subprocess.run(["wget", "-q", f"https://files.rcsb.org/download/{pdb_id}.pdb", "-O", raw_path])
    with open(raw_path) as f:
        lines = f.readlines()
    hetero = set(l[17:20].strip() for l in lines if l.startswith("HETATM"))
    return hetero


def auto_dock_precedent(rule_name, original_smiles, fixed_smiles, use_cache=True):
    """rule_name으로 DOCKING_TARGETS에서 표적 정보를 찾아 도킹 실행.
    pdb_id/ligand_code가 None이면 아직 사람 확인이 안 된 표적이므로
    명확히 에러를 반환한다(자동으로 대충 진행하지 않음)."""
    if rule_name not in DOCKING_TARGETS:
        return {"error": f"'{rule_name}'은 DOCKING_TARGETS에 등록되지 않음"}

    target_info = DOCKING_TARGETS[rule_name]
    if target_info["pdb_id"] is None or target_info["ligand_code"] is None:
        return {"error": f"'{rule_name}' 표적의 pdb_id/ligand_code가 아직 "
                          f"확인 안 됨. inspect_hetatm()으로 먼저 확인 후 "
                          f"DOCKING_TARGETS를 채워주세요."}

    cache = _load_cache() if use_cache else {}
    cache_key = f"{rule_name}|{original_smiles}|{fixed_smiles}"
    if use_cache and cache_key in cache:
        return cache[cache_key]

    if not os.path.exists(VINA_BIN):
        return {"error": f"vina_bin이 {VINA_BIN}에 없음. 다운로드 먼저 진행하세요."}

    receptor_pdbqt, lines = prepare_receptor(target_info["pdb_id"], target_info["keep_hetatm_codes"])

    ligand_lines = [l for l in lines if l.startswith("HETATM") and l[17:20].strip() == target_info["ligand_code"]]
    if not ligand_lines:
        return {"error": f"리간드 코드 {target_info['ligand_code']}를 PDB에서 찾을 수 없음"}
    coords = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in ligand_lines]
    box_center = [sum(c[i] for c in coords) / len(coords) for i in range(3)]

    scores = {}
    for label, smi in [("original", original_smiles), ("fixed", fixed_smiles)]:
        lig_pdbqt = prepare_ligand_pdbqt(smi, f"targets/{rule_name}_{label}")
        if lig_pdbqt is None:
            return {"error": f"{label} 리간드 준비 실패"}
        scores[label] = run_docking_cli(lig_pdbqt, receptor_pdbqt, f"targets/{rule_name}_{label}",
                                          box_center, target_info["box_size"])

    result = {
        "target": target_info["target_name"], "pdb_id": target_info["pdb_id"], "rule": rule_name,
        "score_original": scores["original"], "score_fixed": scores["fixed"],
        "delta": (scores["fixed"] - scores["original"])
                 if scores["original"] is not None and scores["fixed"] is not None else None
    }

    if use_cache:
        cache[cache_key] = result
        _save_cache(cache)

    return result
