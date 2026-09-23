
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

    elif edit_type == "cleave_bond_add_oh":
        # cleave_bond와 달리, 가수분해로 끊어지는 쪽(hetero_idx_in_pattern,
        # 보통 P/S 등 헤테로원자)에는 물의 OH가 새로 결합한다고 보고 산소를
        # 명시적으로 추가한다. 반대쪽(leaving_idx_in_pattern, 보통 O-C의 C측
        # 알코올/페놀이 되는 쪽)은 기존 cleave_bond처럼 암묵적 수소로 채움.
        pair = candidate["cleave_pair_in_pattern"]
        idx1 = match[pair[0]]
        idx2 = match[pair[1]]
        bond = rwmol.GetBondBetweenAtoms(idx1, idx2)
        if bond is None:
            return None
        rwmol.RemoveBond(idx1, idx2)

        hetero_key = candidate.get("hetero_idx_in_pattern")
        hetero_idx = match[hetero_key] if hetero_key is not None else idx1
        leaving_idx = idx2 if hetero_idx == idx1 else idx1

        rwmol.GetAtomWithIdx(leaving_idx).SetNoImplicit(False)

        new_o_idx = rwmol.AddAtom(Chem.Atom(8))
        rwmol.AddBond(hetero_idx, new_o_idx, Chem.BondType.SINGLE)
        rwmol.GetAtomWithIdx(hetero_idx).SetNoImplicit(False)
        new_o_atom = rwmol.GetAtomWithIdx(new_o_idx)
        new_o_atom.SetNoImplicit(False)
        new_o_atom.SetNumExplicitHs(1)

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
        # 실제 RDKit BRENK Aliphatic_long_chain SMARTS는
        # "[R0&D2][R0&D2][R0&D2][R0&D2]" -> 원소 무관, degree==2 연속
        # 4개. degree를 3으로 만드는 분기(branch)로 끊는다.
        # matches[0]만 쓰면 우연히 말단/분기 원자를 앵커로 잡아 사슬
        # 추적이 즉시 끊길 수 있어(n<4), 모든 매치를 후보로 시도해서
        # 실제로 가장 긴 사슬이 나오는 앵커를 채택한다.

        def _is_chain_member(atom):
            return (not atom.GetIsAromatic()) and (not atom.IsInRing()) and atom.GetDegree() == 2

        def _trace_chain(mol, start, avoid):
            chain = []
            current, prev = start, avoid
            while _is_chain_member(mol.GetAtomWithIdx(current)):
                chain.append(current)
                nbs = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors() if n.GetIdx() != prev]
                if len(nbs) != 1:
                    break
                nxt = nbs[0]
                if nxt in chain:
                    break
                prev, current = current, nxt
                if len(chain) > 60:
                    break
            return chain

        chain_pos = candidate["chain_start_idx_in_pattern"]
        chain_atoms = []
        for m in matches:
            cand_start = m[chain_pos]
            cand_atom = mol.GetAtomWithIdx(cand_start)
            cand_neighbors = [n.GetIdx() for n in cand_atom.GetNeighbors()]
            cand_traces = [_trace_chain(mol, nb, cand_start) for nb in cand_neighbors]
            cand_traces.sort(key=len, reverse=True)
            cleft = cand_traces[0] if len(cand_traces) > 0 else []
            cright = cand_traces[1] if len(cand_traces) > 1 else []
            cand_chain = list(reversed(cleft)) + [cand_start] + cright
            if len(cand_chain) > len(chain_atoms):
                chain_atoms = cand_chain

        n = len(chain_atoms)
        if n < 4:
            return None

        insert_positions = list(range(3, n, 4))
        if insert_positions and (n - 1 - insert_positions[-1]) >= 4:
            insert_positions.append(min(n - 1, insert_positions[-1] + 4))
        elif not insert_positions:
            insert_positions = [min(3, n - 1)]

        branch_atomic_num = 6  # 메틸 분기(탄소). 산소로 분기하면 기존 에테르 O 옆에서
                                # O-C-O 패턴이 생겨 'het-C-het_not_in_ring'을 새로 유발함

        def _find_branch_carbon(p):
            for cp in (p, p - 1, p + 1):
                if 0 <= cp < n:
                    a = rwmol.GetAtomWithIdx(chain_atoms[cp])
                    if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0:
                        return chain_atoms[cp]
            for radius in range(2, n):
                for cp in (p - radius, p + radius):
                    if 0 <= cp < n:
                        a = rwmol.GetAtomWithIdx(chain_atoms[cp])
                        if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0:
                            return chain_atoms[cp]
            return None

        added = 0
        used_targets = set()
        for p in insert_positions:
            target_idx = _find_branch_carbon(p)
            if target_idx is None or target_idx in used_targets:
                continue
            new_atom = rwmol.AddAtom(Chem.Atom(branch_atomic_num))
            rwmol.AddBond(target_idx, new_atom, Chem.BondType.SINGLE)
            rwmol.GetAtomWithIdx(target_idx).SetNoImplicit(False)
            used_targets.add(target_idx)
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

        if edit_type not in ("cleave_bond", "cleave_bond_add_oh") and not allow_counterion and '.' in new_smiles:
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
