# Claude Code DeepSeek V4 Pro Speed and Price Results

Method: five independent Claude Code `--print` single-turn sessions using the fixed prompt, no session persistence, and `--tools ""`. The generated programs were not executed. Completion means that the request succeeded and `modelUsage` verified `deepseek-v4-pro[1m]`; format deviations are recorded but do not change completion.

Reference prices: $0.435 per 1M non-cached input tokens and $0.87 per 1M output tokens.

| Run | Model | API time (s) | Input tokens | Output tokens | Output tokens/s | Total cost (USD) | Cost/1k output (USD) | Completed |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | deepseek-v4-pro[1m] | 26.957 | 6,317 | 2,601 | 96.487 | 0.005011 | 0.001926 | Yes |
| 2 | deepseek-v4-pro[1m] | 24.291 | 6,317 | 2,546 | 104.812 | 0.004963 | 0.001949 | Yes |
| 3 | deepseek-v4-pro[1m] | 22.331 | 6,317 | 2,262 | 101.296 | 0.004716 | 0.002085 | Yes |
| 4 | deepseek-v4-pro[1m] | 55.277 | 6,317 | 5,710 | 103.298 | 0.007716 | 0.001351 | Yes |
| 5 | deepseek-v4-pro[1m] | 27.542 | 6,317 | 2,826 | 102.606 | 0.005207 | 0.001842 | Yes |

## Aggregate

| Mean time (s) | Mean input tokens | Mean output tokens | Mean output tokens/s | Speed SD | Mean cost (USD) | Mean cost/1k output (USD) | Completion |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 31.280 | 6,317.0 | 3,189.0 | 101.700 | 3.178 | 0.005522 | 0.001831 | 100% |

Total reference cost for five runs: **$0.027612**.

## Notes

- `api_time_s` is end-to-end Claude Code wall-clock time.
- Model identity was verified as `deepseek-v4-pro[1m]` in every response.
- Runs 1, 3, and 4 used Markdown/prose wrappers.
- Runs 2 and 5 used pseudo `Write` or tool-call wrappers.
- After wrapper removal, all five code bodies passed Python syntax and all required-function checks.

## Code Quality Scores

The uniform 100-point rubric assigns 10 points to raw-output conformance, 15 to required functions and signatures, 35 to numerical agreement with an independent Newmark reference, 20 to valid artifacts and metrics, 15 to assertion tests and direct execution, and 5 to dependency and implementation requirements.

| Run | Score | Main deductions |
|---:|---:|---|
| 1 | 62.0 | Markdown/prose wrapper (-5), unstable and numerically invalid Newmark response (-25), and non-finite CSV/JSON results (-8). |
| 2 | 72.0 | Pseudo `Write` wrapper (-6), materially incorrect Newmark response (-20), and three required metrics absent from `summary.json` (-2). |
| 3 | 64.0 | Markdown fences (-3), unstable and numerically invalid Newmark response (-25), and non-finite CSV/JSON results (-8). |
| 4 | 95.0 | Markdown/prose wrapper instead of code-only output (-5). The extracted code otherwise passed all checks. |
| 5 | 94.0 | Pseudo `Write`/tool-call wrapper instead of code-only output (-6). The extracted code otherwise passed all checks. |

Mean code-quality score: **77.4/100**.
