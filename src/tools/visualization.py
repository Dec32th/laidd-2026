from rdkit import Chem
from rdkit.Chem import Draw, AllChem


def _truncate(text, max_len=60):
    if text is None:
        return ""
    return text if len(text) <= max_len else text[:max_len] + "..."


def visualize_fix_process(loop_result, mols_per_row=3, sub_img_size=(320, 320)):
    """iterative_fix_loop의 결과를 받아, 각 단계의 분자 구조를
    규칙/후보/판단 이유와 함께 2D 그리드 이미지로 반환.
    치환 과정 중 replace_ring(형태 변화가 핵심 논점인 편집)이 쓰였다면,
    해당 단계의 3D 형태 비교(py3Dmol)도 함께 생성해서 반환한다."""
    from src.tools.replacement_library import get_replacement_candidates

    history = loop_result['history']
    mols = []
    legends = []
    shape_relevant_step = None  # replace_ring이 쓰인 첫 스텝을 기록

    for h in history:
        mol = Chem.MolFromSmiles(h['smiles'])
        mols.append(mol)

        if h['step'] == 0:
            legend = "Step 0 (원본)"
        else:
            reason = _truncate(h.get('candidate_reason', ''), 50)
            legend = f"Step {h['step']}: {h['fixed_rule']}\n-> {h['candidate_used']}\n({reason})"

            rule_info = get_replacement_candidates(h['fixed_rule'])
            if rule_info:
                for c in rule_info.get('candidates', []):
                    if c.get('name') == h['candidate_used'] and c.get('edit_type') == 'replace_ring':
                        if shape_relevant_step is None:
                            shape_relevant_step = h['step']

        legends.append(legend)

    img_2d = Draw.MolsToGridImage(
        mols, molsPerRow=mols_per_row, subImgSize=sub_img_size, legends=legends
    )

    result = {"image_2d": img_2d, "shape_relevant_step": shape_relevant_step}

    if shape_relevant_step is not None:
        prev_smiles = history[shape_relevant_step - 1]['smiles']
        curr_smiles = history[shape_relevant_step]['smiles']
        result["shape_comparison_smiles"] = (prev_smiles, curr_smiles)

    return result


def get_3d_mol(smiles, random_seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=random_seed) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def show_3d_shape_comparison(original_smiles, fixed_smiles, style="stick"):
    """치환 전후 두 분자를 3D 인터랙티브 뷰(py3Dmol)로 나란히 비교."""
    import py3Dmol

    mol_o = get_3d_mol(original_smiles)
    mol_f = get_3d_mol(fixed_smiles)
    if mol_o is None or mol_f is None:
        return None

    mb_o = Chem.MolToMolBlock(mol_o)
    mb_f = Chem.MolToMolBlock(mol_f)

    view = py3Dmol.view(width=800, height=400, viewergrid=(1, 2))
    view.addModel(mb_o, 'mol', viewer=(0, 0))
    view.addModel(mb_f, 'mol', viewer=(0, 1))

    if style == "surface":
        view.setStyle({'stick': {}}, viewer=(0, 0))
        view.setStyle({'stick': {}}, viewer=(0, 1))
        view.addSurface(py3Dmol.VDW, {'opacity': 0.6}, viewer=(0, 0))
        view.addSurface(py3Dmol.VDW, {'opacity': 0.6}, viewer=(0, 1))
    else:
        view.setStyle({style: {}}, viewer=(0, 0))
        view.setStyle({style: {}}, viewer=(0, 1))

    view.zoomTo(viewer=(0, 0))
    view.zoomTo(viewer=(0, 1))
    return view


def visualize_fix_process_full(loop_result, save_prefix=None):
    """단계별 2D 구조(사유 포함)를 항상 보여주고, replace_ring처럼 형태
    자체가 핵심인 편집이 포함된 경우 3D 형태 비교까지 함께 보여준다.
    save_prefix가 주어지면 이미지를 파일로 저장(예: 'quinone' ->
    'quinone_2d.png', 'quinone_3d_step2.html')."""
    result = visualize_fix_process(loop_result)

    print("=== 단계별 구조 변화 (2D) ===")
    display(result["image_2d"])

    if save_prefix:
        with open(f"{save_prefix}_2d.png", "wb") as f:
            f.write(result["image_2d"].data)
        print(f"저장됨: {save_prefix}_2d.png")

    if result["shape_relevant_step"] is not None:
        print(f"\n=== Step {result['shape_relevant_step']}: 형태 변화가 핵심인 치환 감지, 3D 비교 ===")
        orig, fixed = result["shape_comparison_smiles"]
        view = show_3d_shape_comparison(orig, fixed, style="surface")
        if view:
            view.show()
            if save_prefix:
                html_str = view._make_html()
                with open(f"{save_prefix}_3d_step{result['shape_relevant_step']}.html", "w") as f:
                    f.write(html_str)
                print(f"저장됨: {save_prefix}_3d_step{result['shape_relevant_step']}.html")

    return result
