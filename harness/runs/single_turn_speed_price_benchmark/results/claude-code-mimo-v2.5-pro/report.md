# Claude Code MiMo V2.5 Pro Speed and Price Results

Method: five independent Claude Code `--print` single-turn sessions using the fixed prompt, no session persistence, and `--tools ""`. The generated programs were not executed. Completion means that the one-turn request succeeded and `modelUsage` verified `mimo-v2.5-pro`; format deviations are recorded but do not change completion.

Reference prices: $0.435 per 1M non-cached input tokens and $0.87 per 1M output tokens. Cached-input details are not recorded.

| Run | Model | API time (s) | Input tokens | Output tokens | Output tokens/s | Total cost (USD) | Cost/1k output (USD) | Completed |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | mimo-v2.5-pro | 53.151 | 6,039 | 2,548 | 47.939 | 0.004844 | 0.001901 | Yes |
| 2 | mimo-v2.5-pro | 50.353 | 1,943 | 2,346 | 46.591 | 0.002886 | 0.001230 | Yes |
| 3 | mimo-v2.5-pro | 47.954 | 1,943 | 2,271 | 47.358 | 0.002821 | 0.001242 | Yes |
| 4 | mimo-v2.5-pro | 56.493 | 1,943 | 2,818 | 49.883 | 0.003297 | 0.001170 | Yes |
| 5 | mimo-v2.5-pro | 60.496 | 1,943 | 3,062 | 50.615 | 0.003509 | 0.001146 | Yes |

## Aggregate

| Mean time (s) | Mean input tokens | Mean output tokens | Mean output tokens/s | Speed SD | Mean cost (USD) | Mean cost/1k output (USD) | Completion |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 53.689 | 2,762.2 | 2,609.0 | 48.477 | 1.706 | 0.003471 | 0.001338 | 100% |

Total reference cost for five runs: **$0.017357**.

## Notes

- `api_time_s` is end-to-end Claude Code wall-clock time.
- Model identity was verified from the response `modelUsage` field.
- The active endpoint during this test was Xiaomi's Token Plan endpoint, so the dollar figures are pay-as-you-go reference estimates rather than subscription charges.
- Four responses used pseudo `Write` tool-call wrappers and one response used Markdown fences. These violations of `Output code only` were retained in the raw evidence and did not change completion under the adopted policy.

## Code Quality Scores

The uniform 100-point rubric assigns 10 points to raw-output conformance, 15 to required functions and signatures, 35 to numerical agreement with an independent Newmark reference, 20 to valid artifacts and metrics, 15 to assertion tests and direct execution, and 5 to dependency and implementation requirements.

| Run | Score | Main deductions |
|---:|---:|---|
| 1 | 94.0 | Pseudo `Write` wrapper instead of code-only output (-6). The extracted code otherwise passed all checks. |
| 2 | 92.8 | Pseudo `Write` wrapper (-6) and missing required `main()` function/interface (-1.2). Numerical and artifact checks passed. |
| 3 | 97.0 | Markdown code fences instead of code-only output (-3). The extracted code otherwise passed all checks. |
| 4 | 74.0 | Pseudo `Write` wrapper (-6) and materially incorrect Newmark response relative to the independent reference (-20). |
| 5 | 94.0 | Pseudo `Write` wrapper instead of code-only output (-6). The extracted code otherwise passed all checks. |

Mean code-quality score: **90.4/100**.
