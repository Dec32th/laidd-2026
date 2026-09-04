# 미팅 예시 2 — Aliphatic_long_chain 완전해결

```
============================================================
치환 감사추적 리포트
============================================================
원본 분자: CCCCCCCCCCCCCCCC
최종 상태: success
최종 분자: CCCCC(C)CCCC(C)CCCC(C)CCC

--- 단계별 이력 ---

[스텝 0] CCCCCCCCCCCCCCCC
  진단된 문제: ['Aliphatic_long_chain']

[스텝 1] CCCCC(C)CCCC(C)CCCC(C)CCC
  고친 규칙: Aliphatic_long_chain (판단 근거: 규칙 기반(충돌 조정: 남는 문제 수 적은 순 - {'Aliphatic_long_chain': 0}))
  적용된 치환: multi-ether chain (multiple O inserted for long chains)
  candidate 선택 근거: 규칙 기반(고정 인덱스 우선, 실패/미해소 시 같은 규칙 내 다른 candidate로 재시도) (candidate_idx=1, 완전 해소)
============================================================
```
