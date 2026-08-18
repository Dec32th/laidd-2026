"""iterative_fix_loop 결과를 사람이 읽기 좋은 감사추적 리포트로 변환.
심사/발표 자료용 — AI가 왜 그렇게 판단했는지, 언제 사람 검토로 넘겼는지를
그대로 보여준다."""


def generate_audit_report(result, original_smiles=None):
    lines = []
    lines.append("=" * 60)
    lines.append("치환 감사추적 리포트")
    lines.append("=" * 60)
    if original_smiles:
        lines.append(f"원본 분자: {original_smiles}")
    lines.append(f"최종 상태: {result['status']}")
    lines.append(f"최종 분자: {result.get('final_smiles', '')}")
    lines.append("")

    lines.append("--- 단계별 이력 ---")
    for h in result.get('history', []):
        step = h.get('step')
        lines.append(f"\n[스텝 {step}] {h.get('smiles', '')}")
        problems = h.get('problems')
        if problems:
            rule_names = [p['rule_name'] for p in problems]
            lines.append(f"  진단된 문제: {rule_names}")
        if 'fixed_rule' in h:
            lines.append(f"  고친 규칙: {h['fixed_rule']} (판단 근거: {h.get('problem_reason', '')})")
            lines.append(f"  적용된 치환: {h.get('candidate_used', '')}")
            lines.append(f"  candidate 선택 근거: {h.get('candidate_reason', '')}")
            debate_attempts = h.get('debate_rounds')
            if debate_attempts:
                lines.append("  --- 토의(debate) 시도 기록 ---")
                for attempt in debate_attempts:
                    lines.append(f"    ▸ {attempt['rule']}[idx={attempt['candidate_idx']}] 최종: {attempt['verdict']}")
                    for r in attempt['rounds']:
                        role = r.get('role')
                        round_num = r.get('round')
                        text = r.get('text', {})
                        if role == 'critic':
                            lines.append(f"      [R{round_num} critic] {text.get('verdict')}: {text.get('reason', '')}")
                        else:
                            lines.append(f"      [R{round_num} proposer] {text.get('stance')}: {text.get('argument', '')}")
    
    if result.get('skipped_rules'):
        lines.append("\n--- 처리 못 하고 넘긴 규칙 (라이브러리 미등록 또는 사람 검토 필요) ---")
        for detail in result.get('skipped_details', []):
            lines.append(f"  - {detail['rule_name']}: {detail['reason']}")

    if result['status'] == 'stuck':
        lines.append(f"\n--- stuck 사유 ---\n{result.get('reason_detail', result.get('reason', ''))}")

    lines.append("=" * 60)
    return "\n".join(lines)
