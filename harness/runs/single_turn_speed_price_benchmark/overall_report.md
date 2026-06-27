# Cross-Model Single-Turn Speed and Price Report

## Summary

Twenty formal requests were completed: five independent single-turn runs for each of four models.

| Model | Mean time (s) | Mean input tokens | Mean output tokens | Mean output tokens/s | Speed SD | Mean cost (USD) | Mean cost/1k output (USD) | Completion |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GPT-5.5 | 57.821 | 10,949.0 | 2,780.6 | 48.129 | 1.327 | 0.138163 | 0.049708 | 100% |
| MiMo V2.5 Pro | 53.689 | 2,762.2 | 2,609.0 | 48.477 | 1.706 | 0.003471 | 0.001338 | 100% |
| MiMo V2.5 Pro UltraSpeed | 6.439 | 3,129.6 | 2,597.4 | 409.751 | 60.107 | 0.010863 | 0.004257 | 100% |
| DeepSeek V4 Pro | 31.280 | 6,317.0 | 3,189.0 | 101.700 | 3.178 | 0.005522 | 0.001831 | 100% |

The total reference cost across all 20 formal requests was **$0.790100**.

## Main Findings

MiMo V2.5 Pro UltraSpeed was the fastest model in this test, averaging 409.751 output tokens/s. It was approximately 4.0 times as fast as DeepSeek V4 Pro and 8.5 times as fast as ordinary MiMo V2.5 Pro or GPT-5.5.

DeepSeek V4 Pro ranked second in speed at 101.700 output tokens/s and showed low run-to-run rate variation despite one longer response.

Ordinary MiMo V2.5 Pro and GPT-5.5 produced nearly identical mean output rates, 48.477 and 48.129 output tokens/s respectively. Their reference costs differed substantially because of their published API prices.

Ordinary MiMo V2.5 Pro had the lowest mean reference cost per 1,000 output tokens, followed by DeepSeek V4 Pro, MiMo V2.5 Pro UltraSpeed, and GPT-5.5.

## Output Conformance

- GPT-5.5 returned code-only outputs in all five runs and passed the static syntax and required-function checks.
- MiMo V2.5 Pro returned four pseudo `Write` wrappers and one Markdown-fenced response. These deviations were recorded without changing the completion flag.
- MiMo V2.5 Pro UltraSpeed returned one pseudo `Write` wrapper and four Markdown-fenced responses. Run 1 misspelled `newmark_average_acceleration` as `newmark_average_accelration`; Runs 2–5 passed the inner-code syntax and required-function checks after wrapper removal.
- DeepSeek V4 Pro returned three Markdown/prose-wrapped responses and two pseudo `Write`/tool-call wrappers. After wrapper removal, all five code bodies passed syntax and required-function checks.

## Interpretation Boundaries

`api_time_s` is end-to-end local CLI wall-clock time. It includes CLI startup, network transfer, provider latency, and response finalization; it is not a server-side time-to-first-token or streaming-only measurement.

The benchmark measures isolated single-turn generation, not agent workflow efficiency. It excludes generated-code execution, file repair, solver execution, iterative debugging, and engineering-result validation.

Input-token accounting is not directly comparable across CLI/provider implementations because their system prompts and automatic prefix caching differ. The speed metric uses each CLI's reported output-token usage divided by measured wall-clock time.

All costs are reference estimates from published standard rates, not invoices. MiMo Token Plan or other subscription accounting may differ.

## Evidence

- [GPT-5.5 report](results/codex-gpt-5.5/report.md)
- [MiMo V2.5 Pro report](results/claude-code-mimo-v2.5-pro/report.md)
- [MiMo V2.5 Pro UltraSpeed report](results/claude-code-mimo-v2.5-pro-ultraspeed/report.md)
- [DeepSeek V4 Pro report](results/claude-code-deepseek-v4-pro/report.md)

Each result directory also contains `results.csv`, `results.json`, `summary.json`, five per-run `result.json` files, raw CLI responses/events, prompts, stderr logs, and generated model output.
