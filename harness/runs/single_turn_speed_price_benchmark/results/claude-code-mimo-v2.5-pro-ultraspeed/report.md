# Claude Code MiMo V2.5 Pro UltraSpeed Speed and Price Results

Method: five independent Claude Code `--print` single-turn sessions using the fixed prompt, no session persistence, and `--tools ""`. The generated programs were not executed. Completion means that the request succeeded and `modelUsage` verified `mimo-v2.5-pro-ultraspeed`; format or content deviations are recorded separately.

Reference prices: $1.305 per 1M non-cached input tokens and $2.61 per 1M output tokens. These are the published limited-time UltraSpeed rates, three times the ordinary MiMo V2.5 Pro rates.

| Run | Model | API time (s) | Input tokens | Output tokens | Output tokens/s | Total cost (USD) | Cost/1k output (USD) | Completed |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | mimo-v2.5-pro-ultraspeed | 5.505 | 6,048 | 2,606 | 473.383 | 0.014694 | 0.005639 | Yes |
| 2 | mimo-v2.5-pro-ultraspeed | 6.544 | 6,048 | 2,278 | 348.089 | 0.013838 | 0.006075 | Yes |
| 3 | mimo-v2.5-pro-ultraspeed | 5.402 | 1,184 | 2,278 | 421.667 | 0.007491 | 0.003288 | Yes |
| 4 | mimo-v2.5-pro-ultraspeed | 6.359 | 1,184 | 2,920 | 459.208 | 0.009166 | 0.003139 | Yes |
| 5 | mimo-v2.5-pro-ultraspeed | 8.386 | 1,184 | 2,905 | 346.410 | 0.009127 | 0.003142 | Yes |

## Aggregate

| Mean time (s) | Mean input tokens | Mean output tokens | Mean output tokens/s | Speed SD | Mean cost (USD) | Mean cost/1k output (USD) | Completion |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6.439 | 3,129.6 | 2,597.4 | 409.751 | 60.107 | 0.010863 | 0.004257 | 100% |

Total reference cost for five runs: **$0.054317**.

## Notes

- `api_time_s` is end-to-end Claude Code wall-clock time.
- Model identity was verified as `mimo-v2.5-pro-ultraspeed` in every response.
- Run 1 used a pseudo `Write` wrapper and misspelled the required function as `newmark_average_accelration`.
- Runs 2–5 used Markdown code fences; Run 5 also included explanatory prose.
- After wrapper removal, Runs 2–5 passed Python syntax and required-function checks.

## Code Quality Scores

The uniform 100-point rubric assigns 10 points to raw-output conformance, 15 to required functions and signatures, 35 to numerical agreement with an independent Newmark reference, 20 to valid artifacts and metrics, 15 to assertion tests and direct execution, and 5 to dependency and implementation requirements.

| Run | Score | Main deductions |
|---:|---:|---|
| 1 | 60.8 | Pseudo `Write` wrapper (-6), misspelled required Newmark function/interface (-1.2), unstable NaN-producing integration (-24), and non-finite CSV/JSON results (-8). |
| 2 | 97.0 | Markdown code fences instead of code-only output (-3). The extracted code otherwise passed all checks. |
| 3 | 67.0 | Markdown fences (-3), unstable NaN-producing Newmark response (-22), and non-finite CSV/JSON results (-8). |
| 4 | 67.0 | Markdown fences (-3), unstable NaN-producing Newmark response (-22), and non-finite CSV/JSON results (-8). |
| 5 | 65.0 | Markdown fences plus explanatory prose (-5), unstable NaN-producing Newmark response (-22), and non-finite CSV/JSON results (-8). |

Mean code-quality score: **71.4/100**.
