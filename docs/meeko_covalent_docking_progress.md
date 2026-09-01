# Meeko + AutoDock-GPU Covalent Docking 진행 기록 (노트북 58)

## 완료된 단계

1. AutoDock-GPU 컴파일: Colab GPU 런타임에서 GPU_INCLUDE_PATH/
   GPU_LIBRARY_PATH를 CUDA 경로로 설정 후 make DEVICE=CUDA 성공.
2. Meeko/ProDy/gemmi 설치 완료.
3. 6JX4 receptor 준비:

   mk_prepare_receptor.py --read_pdb 6jx4_raw.pdb -o 6jx4_receptor -p --delete_residues "A:1101" --allow_bad_res --default_altloc A

   - A:1101(YY3 원본 리간드) 제거 필요 (안 하면 template 매칭 에러)
   - 결정구조 특성상 다수 잔기(697, 708, 709 등)의 곁사슬 원자 결측
     되어 있어 --allow_bad_res로 자동 처리. Cys797은 SG/HG 포함 온전히
     보존됨(좌표 -52.745, -6.415, -17.849, 노트북 54 실측치와 일치).
4. 오시메르티닙 tethered 리간드 준비:
   - SMILES를 RDKit으로 3D 좌표 생성(EmbedMolecule + MMFFOptimizeMolecule)
     후 SDF로 저장 (mk_prepare_ligand.py는 SMILES 텍스트 직접 입력을
     지원하지 않고 SDF/MOL2/MOL만 지원)
   - --receptor 옵션은 PDBQT가 아니라 원본 PDB/mmCIF만 지원함
     (KeyError: 'pdbqt'로 확인)

   mk_prepare_ligand.py -i osimertinib.sdf -o osimertinib_tethered.pdbqt --receptor 6jx4_raw.pdb --rec_residue "A:CYS:797" --tether_smarts "C=CC(=O)N" --tether_smarts_indices 1 2 -v

   - 성공: 출력 PDBQT가 BEGIN_RES CYS A 797 블록(AutoDock 유연잔기 형식)
     으로 시작 — Cys797 사이드체인이 리간드와 함께 하나의 도킹 단위로
     결합된 것으로 보임(정확한 검증은 다음 단계 필요)

## 다음 세션에서 이어갈 것

1. AutoDock-GPU로 이 tethered 리간드 + receptor + 유연잔기 블록을
   어떻게 함께 넘기는지 확인(--flexres 등 관련 플래그 조사)
2. 실제 도킹 실행, 워헤드-SG 거리 확인(5Å 이내면 노트북 54의 원자매핑
   문제를 완전히 우회하고 성공한 것)
3. 성공 확인되면 이타콘산(비특이적 예상)과 비교해서 방법 신뢰도 검증
4. docking.py에 이 경로를 선택적 기능으로 통합할지 결정

## 프로젝트 현황 (2026-09 기준)

JUMP AI 2026 공모전(4th JUMP AI competition, 팀명 MorForge)은 예선
탈락으로 확정됨. 프로젝트는 예선 통과 여부와 무관하게 지속 개선하여
오픈소스 공개 또는 대학원 진학 포트폴리오로 활용할 계획이었으며,
이 방향은 그대로 유지. 9월 10일 경희대 ML 교수님과 미팅 예정.
