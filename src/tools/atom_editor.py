from rdkit import Chem


def apply_atom_edit_from_rule(smiles: str, rule_name: str, candidate_idx: int = 0):
    """replacement_library의 atom_edit 규칙을 이용해 원자/결합/고리 직접 편집을 수행."""
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

    elif edit_type == "reduce_multi_bond":
        pairs = candidate.get("target_pairs_in_pattern", info.get("target_pairs_in_pattern"))
        ring_atoms_pattern = candidate.get("ring_atoms_in_pattern", info.get("ring_atoms_in_pattern"))
        ring_bonds_pattern = candidate.get("ring_bonds_in_pattern", info.get("ring_bonds_in_pattern"))

        for pair in pairs:
            idx_c = match[pair[0]]
            idx_o = match[pair[1]]
            bond = rwmol.GetBondBetweenAtoms(idx_c, idx_o)
            if bond is None:
                return None
            bond.SetBondType(Chem.BondType.SINGLE)
            rwmol.GetAtomWithIdx(idx_o).SetNoImplicit(False)
            rwmol.GetAtomWithIdx(idx_c).SetNumExplicitHs(0)
            rwmol.GetAtomWithIdx(idx_c).SetNoImplicit(False)

        ring_indices = [match[i] for i in ring_atoms_pattern]
        for a in ring_indices:
            rwmol.GetAtomWithIdx(a).SetIsAromatic(True)

        for b1, b2 in ring_bonds_pattern:
            bidx1, bidx2 = match[b1], match[b2]
            rbond = rwmol.GetBondBetweenAtoms(bidx1, bidx2)
            if rbond is None:
                return None
            rbond.SetBondType(Chem.BondType.AROMATIC)
            rbond.SetIsAromatic(True)

    elif edit_type == "replace_multi":
        for sub in candidate["param"]:
            target_idx = match[sub["idx_in_pattern"]]
            atom = rwmol.GetAtomWithIdx(target_idx)
            atom.SetAtomicNum(sub["new_element"])
            atom.SetFormalCharge(sub.get("new_charge", 0))
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)

    elif edit_type == "remove_substituent":
        remove_idx = match[candidate["remove_idx_in_pattern"]]
        upgrade_idx = match[candidate["upgrade_bond_to_idx_in_pattern"]]
        center_idx = match[candidate.get("center_idx_in_pattern", 0)]

        to_remove = set()
        visited = {center_idx}

        stack = [remove_idx]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            to_remove.add(cur)
            for n in mol.GetAtomWithIdx(cur).GetNeighbors():
                if n.GetIdx() not in visited:
                    stack.append(n.GetIdx())

        visited.add(upgrade_idx)
        upgrade_atom = mol.GetAtomWithIdx(upgrade_idx)
        for n in upgrade_atom.GetNeighbors():
            if n.GetIdx() != center_idx and n.GetIdx() not in to_remove:
                stack2 = [n.GetIdx()]
                while stack2:
                    cur2 = stack2.pop()
                    if cur2 in visited:
                        continue
                    visited.add(cur2)
                    to_remove.add(cur2)
                    for n2 in mol.GetAtomWithIdx(cur2).GetNeighbors():
                        if n2.GetIdx() not in visited:
                            stack2.append(n2.GetIdx())

        for ridx in sorted(to_remove, reverse=True):
            rwmol.RemoveAtom(ridx)

        def _adjust3(idx, removed):
            shift = sum(1 for r in removed if r < idx)
            return idx - shift

        center_new = _adjust3(center_idx, to_remove)
        upgrade_new = _adjust3(upgrade_idx, to_remove)

        bond = rwmol.GetBondBetweenAtoms(center_new, upgrade_new)
        if bond is None:
            return None
        bond.SetBondType(Chem.BondType.DOUBLE)
        rwmol.GetAtomWithIdx(center_new).SetNoImplicit(False)
        rwmol.GetAtomWithIdx(upgrade_new).SetNoImplicit(False)

    elif edit_type == "remove_atom":
        remove_idx = match[candidate["remove_idx_in_pattern"]]
        center_idx = match[candidate.get("center_idx_in_pattern", 0)]

        to_remove = set()
        visited = {center_idx}
        stack = [remove_idx]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            to_remove.add(cur)
            for n in mol.GetAtomWithIdx(cur).GetNeighbors():
                if n.GetIdx() not in visited:
                    stack.append(n.GetIdx())

        for ridx in sorted(to_remove, reverse=True):
            rwmol.RemoveAtom(ridx)

        def _adjust4(idx, removed):
            shift = sum(1 for r in removed if r < idx)
            return idx - shift

        center_new = _adjust4(center_idx, to_remove)
        rwmol.GetAtomWithIdx(center_new).SetNoImplicit(False)

    elif edit_type == "cleave_bond":
        pair = candidate["cleave_pair_in_pattern"]
        idx1 = match[pair[0]]
        idx2 = match[pair[1]]
        bond = rwmol.GetBondBetweenAtoms(idx1, idx2)
        if bond is None:
            return None
        rwmol.RemoveBond(idx1, idx2)
        for idx in (idx1, idx2):
            rwmol.GetAtomWithIdx(idx).SetNoImplicit(False)

    elif edit_type == "open_epoxide":
        pair = candidate["break_pair_in_pattern"]
        idx_o = match[pair[0]]
        idx_c_break = match[pair[1]]

        bond = rwmol.GetBondBetweenAtoms(idx_o, idx_c_break)
        if bond is None:
            return None
        rwmol.RemoveBond(idx_o, idx_c_break)

        frag = Chem.MolFromSmiles("O")
        if frag is None:
            return None
        combined = Chem.CombineMols(rwmol.GetMol(), frag)
        rwmol = Chem.RWMol(combined)
        offset = mol.GetNumAtoms()
        rwmol.AddBond(idx_c_break, offset, Chem.BondType.SINGLE)

        rwmol.GetAtomWithIdx(idx_o).SetNoImplicit(False)
        rwmol.GetAtomWithIdx(idx_c_break).SetNoImplicit(False)
        rwmol.GetAtomWithIdx(offset).SetNoImplicit(False)

    elif edit_type == "replace_ring":
        ring_key = candidate.get("ring_atom_indices_in_pattern", info.get("ring_atom_indices_in_pattern"))
        anchor_key = candidate.get("anchor_indices_in_pattern", info.get("anchor_indices_in_pattern"))
        ring_indices = [match[i] for i in ring_key]
        anchor_idx1 = match[anchor_key[0]]
        anchor_idx2 = match[anchor_key[1]]

        # 안전장치: 고리 원자가 anchor 2개 외에 다른 치환기(메틸기 등)를
        # 갖고 있으면, 그 치환기가 고아가 되어 분자가 조각나므로 치환을
        # 거부한다 (다중 BCP 치환 조각화 버그 재발 방지)
        ring_set = set(ring_indices)
        for ridx in ring_indices:
            ratom = mol.GetAtomWithIdx(ridx)
            for n in ratom.GetNeighbors():
                nidx = n.GetIdx()
                if nidx not in ring_set and nidx not in (anchor_idx1, anchor_idx2):
                    return None

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
        # 다중 조각(fragment) 방지: 결과가 여러 개의 분리된 분자로
        # 나뉘었으면 안전장치가 놓친 조각화로 간주해 무효 처리
        if '.' in new_smiles:
            is_valid = False
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
