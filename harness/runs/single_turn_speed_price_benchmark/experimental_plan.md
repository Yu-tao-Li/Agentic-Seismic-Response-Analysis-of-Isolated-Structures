# Single-Turn Output Speed and Price Benchmark Plan

## 1. Objective

This benchmark compares model output speed and reference token cost under the same fixed Python code-generation request.

It is not an end-to-end agent-workflow benchmark. It does not run ABAQUS, OpenSeesPy, MATLAB, or any other external solver. The generated Python programs are preserved but not executed.

## 2. Models

Each model is tested five times:

- `gpt-5.5`, invoked through Codex CLI 0.137.0.
- `mimo-v2.5-pro`, invoked through Claude Code 2.1.112.
- `mimo-v2.5-pro-ultraspeed`, invoked through Claude Code 2.1.112.
- `deepseek-v4-pro[1m]`, invoked through Claude Code 2.1.112.

## 3. Repetition and Isolation

Every formal run must:

- use the same fixed prompt;
- use a new independent single-turn session;
- avoid conversation-context reuse;
- avoid tool execution;
- avoid executing the generated program; and
- preserve the raw model response and reported usage.

Codex runs use `codex exec --ephemeral --json` in a read-only sandbox. Claude Code runs use `claude --print --output-format json --no-session-persistence --tools ""`.

Provider-side automatic prefix caching may still occur. In accordance with the original plan, cached-input token details are not included in the formal result tables, and all recorded input tokens are priced at the ordinary non-cached input rate.

## 4. Fixed Prompt

```text
Generate one complete Python file named `sdof_newmark_demo.py`.

Output code only. Do not use Markdown fences. Do not use external files. Use only Python standard library, NumPy, and Matplotlib.

The file must analyze the following fixed SDOF oscillator:

- Mass: m = 1000.0 kg
- Stiffness: k = 200000.0 N/m
- Damping ratio: zeta = 0.05
- Gravity: g = 9.81 m/s^2
- Time step: dt = 0.01 s
- Duration: T = 20.0 s
- Initial displacement: u0 = 0.0 m
- Initial velocity: v0 = 0.0 m/s

Define the ground acceleration inside the script as:
ag(t) = 0.30*g*sin(2*pi*1.2*t)*exp(-0.15*t)
      + 0.10*g*sin(2*pi*3.0*t)*exp(-0.20*t)

The equation of motion is:
m*u_ddot + c*u_dot + k*u = -m*ag(t)

The script must include:

1. `make_time_vector(dt, duration)`
2. `ground_acceleration(t, g=9.81)`
3. `compute_damping(m, k, zeta)`
4. `newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)`
5. `compute_response_metrics(t, ag, u, v, a_rel, a_abs)`
6. `save_csv(path, t, ag, u, v, a_rel, a_abs)`
7. `save_summary_json(path, metrics)`
8. `plot_response(path, t, ag, u, v, a_abs)`
9. `run_analysis(output_dir="sdof_output")`
10. `main()`

The Newmark method must use beta = 1/4 and gamma = 1/2.

The script must compute and report:
- natural circular frequency
- natural period
- damping coefficient
- peak relative displacement
- peak relative velocity
- peak relative acceleration
- peak absolute acceleration
- peak pseudo base shear k*u
- time of peak displacement

The script must create the output directory and write:
- `response.csv`
- `summary.json`
- `response.png`

Add three simple assertion-based test functions at the bottom:
- `test_zero_ground_motion()`
- `test_output_shapes()`
- `test_damping_positive()`

Do not include TODOs, placeholders, pseudocode, or omitted sections.
The final code should be directly runnable.
```

## 5. Recorded Fields

Each run records:

| Field | Meaning |
|---|---|
| `run_id` | Unique run identifier |
| `model` | Model identifier verified from the invocation or response |
| `api_time_s` | End-to-end CLI wall-clock time in seconds |
| `input_tokens` | Reported non-cached input tokens used by the formal formula |
| `output_tokens` | Reported output tokens |
| `output_tokens_per_s` | `output_tokens / api_time_s` |
| `input_cost` | Reference input-token cost |
| `output_cost` | Reference output-token cost |
| `total_cost` | Sum of reference input and output cost |
| `cost_per_1k_output_tokens` | Reference total cost per 1,000 output tokens |
| `completed` | Successful one-turn response from the verified requested model |
| `notes` | Output-format, static-code, or other anomalies |

## 6. Reference Prices

The prices used were the published standard API rates available at the time of testing:

| Model | Input per 1M tokens | Output per 1M tokens |
|---|---:|---:|
| GPT-5.5 | $5.000 | $30.000 |
| MiMo V2.5 Pro | $0.435 | $0.870 |
| MiMo V2.5 Pro UltraSpeed | $1.305 | $2.610 |
| DeepSeek V4 Pro | $0.435 | $0.870 |

UltraSpeed uses its published limited-time rate of three times the ordinary MiMo V2.5 Pro rate. Token Plan or subscription billing may differ from the reference pay-as-you-go estimates.

## 7. Formulas

```text
output_tokens_per_s = output_tokens / api_time_s
input_cost = input_tokens / 1e6 * input_price_per_million
output_cost = output_tokens / 1e6 * output_price_per_million
total_cost = input_cost + output_cost
cost_per_1k_output_tokens = total_cost / output_tokens * 1000
```

The reported standard deviation is the sample standard deviation across the five per-run output rates.

## 8. Completion and Conformance Policy

At the user's direction, `completed = true` means the request returned successfully in one turn and the requested model identity was verified. Markdown fences, pseudo tool-call wrappers, prose surrounding the code, or static code defects are disclosed in `notes` and the reports but do not change the completion flag.

This definition should not be confused with strict prompt conformance or proven runtime correctness.
