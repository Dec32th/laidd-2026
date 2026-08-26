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
        "ligand_code": "DNC", "keep_hetatm_codes": ["MG", "SAM"],
        "box_size": [20, 20, 20],
        "caveat": None,  # 매치된 분자가 실제로 카테콜 골격을 가지므로 표적 특이성 근거 있음
    },
    "Michael_acceptor_1": {
        "target_name": "EGFR", "pdb_id": "6JX4",
        "ligand_code": "YY3", "keep_hetatm_codes": [],
        "box_size": [20, 20, 20],
        "caveat": ("[방법론적 한계] 이 도킹은 AutoDock Vina 표준 도킹으로, 공유결합 형성을 전혀 "
                   "모델링하지 못합니다 — 워헤드가 실제 반응 부위(예: 시스테인 잔기)와 형성하는 "
                   "공유결합의 결합력은 반영되지 않고, 오직 비공유 상호작용만 점수화됩니다. "
                   "따라서 warhead 제거 전후 도킹 점수 차이가 작다는 것이(예: <1 kcal/mol) "
                   "'약효 손실이 없다'는 증거가 될 수 없습니다 — 애초에 이 도킹 방법으로는 "
                   "공유결합 억제 메커니즘 자체를 평가할 수 없기 때문입니다. 추가로, EGFR은 "
                   "Michael_acceptor_1 규칙 자체와 특이적 연관이 없는 벤치마크 표적이라, 이 분자가 "
                   "실제로 EGFR을 겨냥한다는 근거도 없습니다. 이 도킹 수치는 참고 이상의 근거로 "
                   "쓰지 말고, warhead 제거 여부는 이 분자가 실제 공유결합 표적을 가진다는 "
                   "별도 증거(precedent, 알려진 반응 부위 등)를 기준으로 판단하세요."),
    },
    "hydroquinone": {
        "target_name": "NQO1", "pdb_id": "1DXO",
        "ligand_code": "DQN", "keep_hetatm_codes": ["FAD"],
        "box_size": [20, 20, 20],
        "caveat": None,
    },
    "quinone_A(370)": {
        "target_name": "NQO1", "pdb_id": "1DXO",
        "ligand_code": "DQN", "keep_hetatm_codes": ["FAD"],
        "box_size": [20, 20, 20],
        "caveat": None,
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


def prepare_ligand_pdbqt(smiles, filename, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
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
        cached_result = dict(cache[cache_key])
        cached_result["caveat"] = target_info.get("caveat")  # caveat은 캐시 무시, 항상 최신 반영
        return cached_result

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
                 if scores["original"] is not None and scores["fixed"] is not None else None,
        "caveat": target_info.get("caveat"),
    }

    if use_cache:
        cache[cache_key] = result
        _save_cache(cache)

    return result

# Michael acceptor류 반응성/선택성에 대한 문헌 기반 참고표.
# 좌표 계산 없이, 알려진 화학 지식을 critic 판단에 참고자료로 제공하기 위함
# (노트북 54: 거리 기반 기하학적 근사는 원자매핑 신뢰성 문제로 폐기, 저비용 대안).
WARHEAD_REACTIVITY_REFERENCE = {
    "Michael_acceptor_1": {
        "warhead_class": "α,β-불포화 카르보닐 (Michael acceptor)",
        "reactivity_note": (
            "이 워헤드 계열의 반응성은 구조에 따라 크게 갈립니다. (1) 아크릴아마이드형 "
            "(예: 오시메르티닙, 아파티닙, 이브루티닙의 warhead)은 특정 시스테인(EGFR Cys797, "
            "BTK Cys481 등)과 선택적으로 반응하도록 설계된 잘 확립된 공유결합 억제제 warhead입니다. "
            "(2) 반면 이타콘산(itaconate) 같은 단순 α,β-불포화 카르복실산도 '반응성이 없다'고 "
            "단정할 수 없습니다 — 이타콘산은 실제로 생체 내에서 KEAP1 Cys151, GAPDH Cys22 등과 "
            "Michael 부가반응으로 공유결합하는 내인성 항염증 대사물질로 잘 알려져 있습니다 "
            "(Mills et al., Nature 2018). 즉 '단순해 보이는 Michael acceptor니까 안전하게 "
            "제거해도 된다'는 가정 자체가 항상 성립하지 않습니다."
        ),
        "practical_guidance": (
            "warhead 제거를 판단할 때는, 도킹 델타가 작다는 것보다 (a) 이 분자가 알려진 "
            "약물/대사물질과 구조적으로 유사한지, (b) 선례 라이브러리에 해당 사례가 있는지를 "
            "우선 근거로 삼으세요. 확실한 근거가 없으면 사람 검토로 넘기는 것이 안전합니다."
        ),
    },
}


def get_warhead_reference(rule_name):
    """규칙 이름으로 반응성 참고표 항목을 조회. 없으면 None."""
    return WARHEAD_REACTIVITY_REFERENCE.get(rule_name)

WARHEAD_REACTIVITY_REFERENCE["alkyl_halide"] = {
    "warhead_class": "알킬 할라이드 (친전자성 알킬화제)",
    "reactivity_note": (
        "알킬 할라이드는 SN2 메커니즘으로 시스테인 등 친핵체와 공유결합을 형성하는 "
        "잘 알려진 친전자성 워헤드입니다. 클로로아세트아미드(chloroacetamide) 계열은 "
        "표적 공유결합 억제제(targeted covalent inhibitor) 설계에 흔히 쓰이는 반응기이며, "
        "질소 머스타드(chlorambucil, cyclophosphamide 등)는 이 반응성을 이용해 DNA를 "
        "알킬화하는 항암제로 설계된 사례입니다. 즉 반응성 자체가 결함이 아니라 의도된 "
        "약효 메커니즘인 경우가 있습니다."
    ),
    "practical_guidance": (
        "할라이드를 단순 제거하기 전에, 이 분자가 알려진 알킬화제 약물/항암제 계열과 "
        "구조적으로 유사한지 먼저 확인하세요. 유사성이 있으면 반응성 제거가 약효 상실로 "
        "이어질 수 있으므로 사람 검토를 권장합니다."
    ),
}

WARHEAD_REACTIVITY_REFERENCE["disulphide"] = {
    "warhead_class": "이황화결합 (다이설파이드, 티올-다이설파이드 교환반응)",
    "reactivity_note": (
        "이황화결합은 티올-다이설파이드 교환반응을 통해 단백질 시스테인과 동적 공유결합을 "
        "형성할 수 있습니다. 오라노핀(auranofin, 항류마티스제)은 티오레독신 환원효소와의 "
        "이 교환반응을 이용해 작용하는 승인 약물이며, 디설피람(disulfiram)도 알데하이드"
        "탈수소효소와 유사한 반응성으로 작용합니다. 즉 이 결합의 반응성이 오히려 약효의 "
        "핵심 메커니즘인 경우가 있습니다."
    ),
    "practical_guidance": (
        "이황화결합을 단순 절단하기 전에, 이 분자가 오라노핀/디설피람류처럼 그 반응성 "
        "자체를 이용하는 약물 계열과 유사한지 확인하세요. 확실한 근거가 없으면 사람 검토로 "
        "넘기는 것이 안전합니다."
    ),
}
