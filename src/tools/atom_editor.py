
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
                continue
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
        rwmol.GetAtomWithIdx(center_new).SetFormalCharge(0)

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

    elif edit_type == "insert_atom":
        pair = candidate["insert_pair_in_pattern"]
        idx1 = match[pair[0]]
        idx2 = match[pair[1]]

        # 과산화물(O-O) 생성 방지: 삽입 위치 양쪽이 이미 산소면 거부
        if rwmol.GetAtomWithIdx(idx1).GetSymbol() == 'O' or rwmol.GetAtomWithIdx(idx2).GetSymbol() == 'O':
            return None

        bond = rwmol.GetBondBetweenAtoms(idx1, idx2)
        if bond is None:
            return None
        rwmol.RemoveBond(idx1, idx2)

        new_atom = Chem.Atom(candidate["param"])
        new_idx = rwmol.AddAtom(new_atom)
        rwmol.AddBond(idx1, new_idx, Chem.BondType.SINGLE)
        rwmol.AddBond(new_idx, idx2, Chem.BondType.SINGLE)

        rwmol.GetAtomWithIdx(idx1).SetNoImplicit(False)
        rwmol.GetAtomWithIdx(idx2).SetNoImplicit(False)

    elif edit_type == "insert_atom_multi_chain":
        # 긴 지방족 사슬 전용(탄소 또는 비카르보닐 에테르 산소로 구성된
        # 사슬 모두 인식): 매치 시작점에서 양쪽 방향을 모두 추적해 더 긴
        # 쪽을 진짜 사슬로 채택한 뒤, 필요한 만큼 산소를 균등 삽입
        start_idx = match[candidate["chain_start_idx_in_pattern"]]

        def _is_chain_member(atom):
            if atom.GetSymbol() == 'C' and not atom.GetIsAromatic():
                return True
            if atom.GetSymbol() == 'O' and atom.GetDegree() == 2 and not atom.GetIsAromatic():
                for nb in atom.GetNeighbors():
                    for bond in nb.GetBonds():
                        if bond.GetBondTypeAsDouble() == 2.0 and nb.GetSymbol() == 'C':
                            other = bond.GetOtherAtom(nb)
                            if other.GetSymbol() == 'O':
                                return False
                return True
            return False

        def _trace_chain(mol, start, avoid):
            chain = [start]
            current = start
            prev = avoid
            while True:
                atom_cur = mol.GetAtomWithIdx(current)
                if not _is_chain_member(atom_cur):
                    break
                next_candidates = [n.GetIdx() for n in atom_cur.GetNeighbors()
                                    if n.GetIdx() != prev and _is_chain_member(n)]
                if not next_candidates:
                    break
                prev, current = current, next_candidates[0]
                chain.append(current)
                if len(chain) > 30:
                    break
            return chain

        start_atom = mol.GetAtomWithIdx(start_idx)
        neighbor_options = [n.GetIdx() for n in start_atom.GetNeighbors() if _is_chain_member(n)]

        best_chain = [start_idx]
        for nb in neighbor_options:
            candidate_chain = [start_idx] + _trace_chain(mol, nb, start_idx)
            if len(candidate_chain) > len(best_chain):
                best_chain = candidate_chain

        chain_atoms = best_chain
        if len(chain_atoms) < 4:
            return None

        n = len(chain_atoms)
        num_inserts = max(1, (n - 1) // 3)
        step = n / (num_inserts + 1)
        insert_positions = sorted(set(int(round(step * (i + 1))) for i in range(num_inserts)))
        insert_positions = [p for p in insert_positions if 0 < p < n]

        insert_after = [chain_atoms[p - 1] for p in insert_positions]
        if not insert_after:
            return None

        added = 0
        for a_idx in insert_after:
            a_pos = chain_atoms.index(a_idx)
            b_idx = chain_atoms[a_pos + 1]
            # 삽입 지점 양쪽이 이미 산소면(과산화물 방지) 건너뜀
            if mol.GetAtomWithIdx(a_idx).GetSymbol() == 'O' or mol.GetAtomWithIdx(b_idx).GetSymbol() == 'O':
                continue
            bond = rwmol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None:
                continue
            rwmol.RemoveBond(a_idx, b_idx)
            new_o = rwmol.AddAtom(Chem.Atom(8))
            rwmol.AddBond(a_idx, new_o, Chem.BondType.SINGLE)
            rwmol.AddBond(new_o, b_idx, Chem.BondType.SINGLE)
            rwmol.GetAtomWithIdx(a_idx).SetNoImplicit(False)
            rwmol.GetAtomWithIdx(b_idx).SetNoImplicit(False)
            added += 1

        if added == 0:
            return None

    elif edit_type == "replace_ring":
        ring_key = candidate.get("ring_atom_indices_in_pattern", info.get("ring_atom_indices_in_pattern"))
        anchor_key = candidate.get("anchor_indices_in_pattern", info.get("anchor_indices_in_pattern"))
        ring_indices = [match[i] for i in ring_key]
        anchor_idx1 = match[anchor_key[0]]
        anchor_idx2 = match[anchor_key[1]]

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
        allow_counterion = candidate.get("allow_counterion", False)
        allow_aromatic_zero_h = candidate.get("allow_aromatic_zero_h", False)

        if edit_type != "cleave_bond" and not allow_counterion and '.' in new_smiles:
            is_valid = False
        for atom in check_mol.GetAtoms():
            if allow_aromatic_zero_h and atom.GetIsAromatic() and atom.GetSymbol() == 'N':
                continue
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
