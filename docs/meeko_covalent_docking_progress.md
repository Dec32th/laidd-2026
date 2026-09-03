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

## 결과: Tethered Docking 검증 완료 (성공)

전체 파이프라인(receptor 준비 → autogrid4 → AutoDock-GPU tethered
docking)이 정상 작동함을 확인. 오시메르티닙(실제 EGFR Cys797 공유결합
억제제) vs 이타콘산(공유결합 근거 없는 대조군)으로 비교:

| 분자 | Best Inter+Intra Energy (kcal/mol) |
|---|---|
| 오시메르티닙 | -5.31 (안정) |
| 이타콘산 | +7.06 (매우 불안정) |

12.4 kcal/mol의 뚜렷한 차이로 방법의 변별력 확인. 노트북 54에서 실패한
원자 매핑 문제를 Meeko의 명시적 원자 이름 태깅(REMARK SMILES IDX)으로
완전히 우회함 — 인덱스 추측이 아니라 Meeko가 보장하는 "결합차수 손실
없는" 왕복 변환을 그대로 신뢰.

## 알아낸 주요 세부사항 (재현을 위한 기록)

- receptor 준비 시 원본 리간드(YY3 등) 반드시 --delete_residues로 제거
- 결정구조의 결측 곁사슬은 --allow_bad_res로 자동 처리(Cys797 자체는
  온전했음)
- --box_center는 원본 리간드 좌표 평균으로 직접 계산 필요
  (--box_enveloping 등 대안 있으나 이번엔 수동 계산 사용)
- GPF 자동생성 시 Si/B 원자타입 관련 파싱 버그 발견 — ligand_types와
  map 줄을 15개로 정확히 세었는데도 autogrid4가 "14개raa"로 오인식.
  Si/B 원자가 필요 없는 경우 해당 줄 제거로 우회
- autogrid4는 apt 별도 패키지(`autogrid`, `autodock`과 분리됨) 설치 필요
- --lfile에 완전히 빈 파일은 허용 안 됨 — 최소 1원자 더미 PDBQT 필요
  (기존 유효 PDBQT에서 원자 1줄의 정확한 컬럼 포맷을 복사해 사용하는
  것이 안전)
- 도킹 결과(DLG)에서 실제 포즈 파싱 시 Meeko의 PDBQTMolecule.from_file
  (is_dlg=True)이 이 AutoDock-GPU 버전 출력과 호환 안 됨 — DOCKED:
  접두어를 직접 제거해 표준 PDBQT로 재구성하는 방식으로 우회
- Tethered 리간드 PDBQT의 BEGIN_RES 블록은 CA/CB까지만 포함하고 SG는
  없음 — "거리 계산"이 아니라 "이 결합형태의 에너지"로 검증해야 함

## 다음 단계

1. docking.py에 선택적 기능(예: `try_covalent_tethered_check`)으로
   통합할지 결정 — 이번엔 Colab 환경 설정(AutoDock-GPU 컴파일,
   autogrid 설치)이 매 세션 필요해 운영 비용이 상당함을 감안
2. Michael_acceptor_1 외 다른 워헤드 규칙에도 같은 방식 적용 가능성
   검토
3. 통합하지 않더라도, 이 성공 사례 자체를 교수님 미팅에서 "정직한
   실패 후 재도전으로 해결한 사례"로 제시 가능
"""

!cd /content/laidd-2026 && git add docs/meeko_covalent_docking_progress.md && git commit -m "Successfully validate Meeko+AutoDock-GPU tethered covalent docking: osimertinib (-5.31 kcal/mol) vs itaconic acid (+7.06 kcal/mol) shows clear discrimination" && git push
