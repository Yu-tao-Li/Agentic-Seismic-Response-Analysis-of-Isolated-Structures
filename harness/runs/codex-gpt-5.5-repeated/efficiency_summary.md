# Codex GPT-5.5 Repeated Run Efficiency Summary

This repeated experiment used the same harness-plus-BRAIN-prompt format for Examples 1~6 and generated fresh outputs under `harness/runs/codex-gpt-5.5-repeated`. No generated outputs from the original `codex-gpt-5.5` run or other agents were copied.

`num_agent_turns` is 1 for each Codex CLI invocation because `codex exec` is a one-shot non-interactive run. For efficiency comparison, `codex_efficiency_action_count` is more useful: it counts completed Codex CLI item events. Example 6 is recorded as an initial failed repeated run followed by one neutral tool-location reminder; the final generated outputs then reached strict three-software agreement.

| Example | Final result | Total duration (min) | Initial duration (min) | Follow-up duration (min) | Codex exec invocations | Completed events | Completed command executions | Agent messages | Manual feedback cycles |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | PASS | 15.812 | 15.812 | 0 | 1 | 96 | 62 | 29 | 0 |
| 2 | PASS | 16.03 | 16.03 | 0 | 1 | 103 | 62 | 35 | 0 |
| 3 | PASS | 25.346 | 25.346 | 0 | 1 | 119 | 68 | 39 | 0 |
| 4 | PASS | 30.135 | 30.135 | 0 | 1 | 130 | 76 | 44 | 0 |
| 5 | PASS | 25.489 | 25.489 | 0 | 1 | 148 | 85 | 41 | 0 |
| 6 | PASS | 124.276 | 25.586 | 98.69 | 2 | 513 | 277 | 166 | 1 |

Examples 1~5 achieved strict three-software agreement in a single repeated Codex invocation. Example 6 first failed because the ABAQUS user-subroutine path was treated as unavailable. A neutral reminder was then issued with the local tool locations, and the follow-up run compiled the ABAQUS UEL, generated ABAQUS nonlinear response output, and updated `cross_validation_report.json` to `strict_three_software_pass: true`.

Example 6 follow-up records:

- `example6/brain/agent_io/codex_followup_tool_hint_20260625_prompt.txt`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_stdout.jsonl`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_stderr.txt`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_last_message.txt`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_process.json`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_tool_locations.txt`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_oneapi_check.log`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_vs_check.log`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_result.json`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_events.json`
- `example6/brain/agent_io/codex_followup_tool_hint_20260625_summary.md`

Machine-readable records:

- `efficiency_summary.csv`
- `efficiency_summary.json`
- `example*/brain/run_metadata.json`
- `example*/brain/agent_io/codex_input.txt`
- `example*/brain/agent_io/codex_stdout.jsonl`
- `example*/brain/agent_io/codex_events.json`
- `example*/brain/agent_io/codex_result.json`
- `example*/brain/cross_validation_report.json`
- `example*/brain/report.md`
