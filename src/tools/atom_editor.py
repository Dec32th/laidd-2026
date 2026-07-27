from rdkit import Chem


def apply_atom_edit_from_rule(smiles: str, rule_name: str, candidate_idx: int = 0):
    """replacement_library의 atom_edit 규칙을 이용해 원자/결합/고리 직접 편집을 수행.
    candidate마다 다른 edit_type을 가질 수 있음 (예: 같은 문제에 대해
    작은 변화(치환기 하나 추가)와 큰 변화(고리 전체 교체)를 후보로 병렬 제시)."""
    from src.tools.replacement_library import get_replacement_candidates
    info = get_replacement_candidates(rule_name)
    if info is None or info.get("edit_method") != "atom_edit":
        return None
    if candidate_idx >= len(info["candidates"]):
        return None

    candidate = info["candidates"][candidate_idx]
    smarts = info["problem_smarts"]

    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts(smarts)
    if mol is None or pattern is None:
        return None

    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return None
    match = matches[0]

    rwmol = Chem.RWMol(mol)
    edit_type = candidate["edit_type"]

    if edit_type == "replace_element":
        target_idx = match[candidate.get("target_idx_in_pattern", info.get("target_idx_in_pattern"))]
        atom = rwmol.GetAtomWithIdx(target_idx)
        atom.SetAtomicNum(candidate["param"])

    elif edit_type == "add_substituent":
        target_idx = match[candidate.get("target_idx_in_pattern", info.get("target_idx_in_pattern"))]
        frag = Chem.MolFromSmiles(candidate["param"])
        if frag is None:
            return None
        combined = Chem.CombineMols(rwmol.GetMol(), frag)
        rwmol = Chem.RWMol(combined)
        offset = mol.GetNumAtoms()
        rwmol.AddBond(target_idx, offset, Chem.BondType.SINGLE)
        atom = rwmol.GetAtomWithIdx(target_idx)
        if atom.GetNumExplicitHs() > 0:
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
        else:
            atom.SetNoImplicit(False)

    elif edit_type == "reduce_bond":
        pair = candidate.get("target_idx_pair_in_pattern", info.get("target_idx_pair_in_pattern"))
        idx1 = match[pair[0]]
        idx2 = match[pair[1]]
        bond = rwmol.GetBondBetweenAtoms(idx1, idx2)
        if bond is None:
            return None
        bond.SetBondType(Chem.BondType.SINGLE)
        for idx in (idx1, idx2):
            atom = rwmol.GetAtomWithIdx(idx)
            atom.SetNoImplicit(False)

    elif edit_type == "replace_multi":
        for sub in candidate["param"]:
            target_idx = match[sub["idx_in_pattern"]]
            atom = rwmol.GetAtomWithIdx(target_idx)
            atom.SetAtomicNum(sub["new_element"])
            atom.SetFormalCharge(sub.get("new_charge", 0))
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)

    elif edit_type == "replace_ring":
        ring_key = candidate.get("ring_atom_indices_in_pattern", info.get("ring_atom_indices_in_pattern"))
        anchor_key = candidate.get("anchor_indices_in_pattern", info.get("anchor_indices_in_pattern"))
        ring_indices = [match[i] for i in ring_key]
        anchor_idx1 = match[anchor_key[0]]
        anchor_idx2 = match[anchor_key[1]]

        anchor1_ring_neighbor = None
        anchor2_ring_neighbor = None
        for ridx in ring_indices:
            ratom = mol.GetAtomWithIdx(ridx)
            neighbor_idxs = [n.GetIdx() for n in ratom.GetNeighbors()]
            if anchor_idx1 in neighbor_idxs:
                anchor1_ring_neighbor = ridx
            if anchor_idx2 in neighbor_idxs:
                anchor2_ring_neighbor = ridx

        if anchor1_ring_neighbor is None or anchor2_ring_neighbor is None:
            return None

        frag = Chem.MolFromSmiles(candidate["param"])
        if frag is None:
            return None

        for ridx in sorted(ring_indices, reverse=True):
            rwmol.RemoveAtom(ridx)

        def _adjust(idx, removed):
            shift = sum(1 for r in removed if r < idx)
            return idx - shift

        anchor_idx1_new = _adjust(anchor_idx1, ring_indices)
        anchor_idx2_new = _adjust(anchor_idx2, ring_indices)

        combined = Chem.CombineMols(rwmol.GetMol(), frag)
        rwmol2 = Chem.RWMol(combined)
        offset = rwmol.GetMol().GetNumAtoms()

        frag_attach1 = None
        frag_attach2 = None
        for atom in frag.GetAtoms():
            if atom.GetSymbol() == '*':
                map_num = atom.GetAtomMapNum()
                if map_num == 1:
                    frag_attach1 = atom.GetIdx() + offset
                elif map_num == 2:
                    frag_attach2 = atom.GetIdx() + offset

        if frag_attach1 is None or frag_attach2 is None:
            return None

        dummy1 = rwmol2.GetAtomWithIdx(frag_attach1)
        dummy2 = rwmol2.GetAtomWithIdx(frag_attach2)
        real_neighbor1 = dummy1.GetNeighbors()[0].GetIdx()
        real_neighbor2 = dummy2.GetNeighbors()[0].GetIdx()

        rwmol2.AddBond(anchor_idx1_new, real_neighbor1, Chem.BondType.SINGLE)
        rwmol2.AddBond(anchor_idx2_new, real_neighbor2, Chem.BondType.SINGLE)
        rwmol2.RemoveAtom(max(frag_attach1, frag_attach2))
        rwmol2.RemoveAtom(min(frag_attach1, frag_attach2))

        rwmol = rwmol2
    else:
        return None

    try:
        new_mol = rwmol.GetMol()
        Chem.SanitizeMol(new_mol)
    except Exception:
        return None

    new_smiles = Chem.MolToSmiles(new_mol)

    check_mol = Chem.MolFromSmiles(new_smiles)
    is_valid = check_mol is not None
    if is_valid:
        for atom in check_mol.GetAtoms():
            if (atom.GetNoImplicit() and atom.GetFormalCharge() == 0
                    and atom.GetSymbol() in ('C', 'N', 'O')
                    and atom.GetTotalNumHs() == 0 and atom.GetDegree() < 4):
                is_valid = False
                break

    return {
        "new_smiles": new_smiles,
        "candidate_used": candidate["name"],
        "rationale": candidate["rationale"],
        "is_valid": is_valid,
    }
