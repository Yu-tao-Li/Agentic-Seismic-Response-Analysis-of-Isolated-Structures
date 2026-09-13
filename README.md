# Harness Evidence Package

This folder contains the program and evidence files used to support the
manuscript conclusions on coding-agent harness evaluation for seismic
analysis programming of isolated structures.

Path-related strings in the text files have been sanitized or converted to
package-relative paths for public release.

## Included Contents

- `harness/prompts/`: BRAIN task prompts for Examples 1-6.
- `harness/runs/codex-gpt-5.5/`: final Codex GPT-5.5 evidence for Examples 1-6.
- `harness/runs/claude-code-deepseek-v4-pro/`: final Claude Code + DeepSeek V4 Pro evidence for Examples 1-6.
- `harness/runs/codex-gpt-5.5-repeated/`: repeated Codex GPT-5.5 run supporting the repeatability discussion.
- `harness/runs/claude-code-mimo-v2.5-pro/example1/` and `harness/runs/claude-code-mimo-v2.5-pro-ultraspeed/example1/`: Example 1 workflow records used for the model-speed comparison.
- `harness/runs/single_turn_speed_price_benchmark/`: controlled five-run single-turn speed, cost, and correctness benchmark.
- `harness/runs/codex-gpt-5.5-repeated/example5/brain/reviewer1_round2_validation/`: canonical corrected Example 5 solver, dimensionless tolerance criteria, 27-run sensitivity analysis, and comparison with the independent Cui et al. reference output prepared for the second revision.
- `external_reference_validation/`: consolidated checks using the archived author-confirmed outputs from Chen's compiled program for Examples 1 and 2, the archived Cui et al. outputs for Example 4, and both the retained historical records and the independently executed published-algorithm output for Example 5. The revised Example 5 series and the external-reference/harness SHA-256 mapping are deposited alongside the comparison evidence.

The retained run folders include the generated programs, input files,
machine-readable outputs, cross-validation reports, run-level reports, plotting
scripts, and ABAQUS raw solver artifacts such as `.odb`, `.msg`, `.prt`, `.dat`,
`.com`, and `.sta` where they directly support the reported evidence chain.

## Notes

MATLAB, OpenSeesPy, and ABAQUS are required to rerun the full structural-analysis
workflows. The stored outputs are included so the reported comparisons can be
audited without rerunning every solver job. The source material available for
the Cui et al. implementations does not establish individual code authorship,
so this package does not attribute those program files to a specific coauthor.
The copyrighted book listing and page scans are not redistributed; the public
Example 5 addition contains numerical response histories, provenance metadata,
and a verification script only.
