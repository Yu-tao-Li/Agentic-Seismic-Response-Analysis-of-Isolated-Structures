# Codex GPT-5.5 Speed and Price Results

Method: five independent `codex exec --ephemeral --json` single-turn sessions using the fixed prompt. The generated programs were not executed. Completion required valid usage, no tool calls, valid Python syntax, and all required functions.

Reference prices: $5.00 per 1M input tokens and $30.00 per 1M output tokens. Cached-input details are not recorded, and all reported input tokens are priced at the standard input rate.

| Run | Model | API time (s) | Input tokens | Output tokens | Output tokens/s | Total cost (USD) | Cost/1k output (USD) | Completed |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | gpt-5.5 | 55.408 | 10,949 | 2,719 | 49.072 | 0.136315 | 0.050134 | Yes |
| 2 | gpt-5.5 | 62.877 | 10,949 | 2,906 | 46.218 | 0.141925 | 0.048839 | Yes |
| 3 | gpt-5.5 | 58.029 | 10,949 | 2,861 | 49.303 | 0.140575 | 0.049135 | Yes |
| 4 | gpt-5.5 | 56.368 | 10,949 | 2,749 | 48.769 | 0.137215 | 0.049915 | Yes |
| 5 | gpt-5.5 | 56.424 | 10,949 | 2,668 | 47.285 | 0.134785 | 0.050519 | Yes |

## Aggregate

| Mean time (s) | Mean input tokens | Mean output tokens | Mean output tokens/s | Speed SD | Mean cost (USD) | Mean cost/1k output (USD) | Completion |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 57.821 | 10,949.0 | 2,780.6 | 48.129 | 1.327 | 0.138163 | 0.049708 | 100% |

Total reference cost for five runs: **$0.690815**.

## Notes

- `api_time_s` is end-to-end CLI wall-clock time.
- `output_tokens` comes from `turn.completed.usage.output_tokens` and includes reported reasoning output tokens.
- Codex adds its own system instructions, so input usage exceeds the fixed user prompt alone.
- All five outputs were code-only, used no tools, parsed as Python, and included all required functions.

## Code Quality Scores

The uniform 100-point rubric assigns 10 points to raw-output conformance, 15 to required functions and signatures, 35 to numerical agreement with an independent Newmark reference, 20 to valid artifacts and metrics, 15 to assertion tests and direct execution, and 5 to dependency and implementation requirements.

| Run | Score | Main deductions |
|---:|---:|---|
| 1 | 100.0 | None. Code-only output, reference-level numerical agreement, valid artifacts, and all tests passed. |
| 2 | 100.0 | None. Code-only output, reference-level numerical agreement, valid artifacts, and all tests passed. |
| 3 | 100.0 | None. Code-only output, reference-level numerical agreement, valid artifacts, and all tests passed. |
| 4 | 100.0 | None. Code-only output, reference-level numerical agreement, valid artifacts, and all tests passed. |
| 5 | 100.0 | None. Code-only output, reference-level numerical agreement, valid artifacts, and all tests passed. |

Mean code-quality score: **100.0/100**.
