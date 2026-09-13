# External Reference Validation

This directory consolidates the external-reference checks for Examples 1, 2, 4, and 5. The archived series for Examples 1, 2, and 4, together with the previously deposited Example 5 series, remain immutable comparison and integrity records. A separate Example 5 reference series was obtained by retranscribing and independently executing the published Cui et al. algorithm in MATLAB R2025b; its numerical outputs are deposited under `example5_cui/independent_published_algorithm_output/`. The cited Cui source is the coauthored book, and no individual code authorship is assigned here.

## Coverage

| Example | External reference | Coverage |
| --- | --- | --- |
| 1 | Chen compiled MDOF program, article `id=408` | Frequencies, periods, and sign-aligned mode shapes |
| 2 | Chen compiled MDOF program, article `id=408` | Top-displacement history |
| 3 | Not available | Internal three-software consistency only |
| 4 | Cui et al. published MATLAB implementation | Top displacement, top relative acceleration, and isolation displacement |
| 5 | Cui et al. published MATLAB implementation | Top response, isolation displacement and restoring force, and hysteresis-loop area |
| 6 | Not available | Internal three-software consistency only |

## Reproduction

Run the following command from this directory with the verified `mypy` environment:

```powershell
D:\anaconda3\envs\mypy\python.exe .\compare_external_references.py
```

The command recomputes the archived comparison metrics from the deposited files and regenerates:

- `external_reference_summary.csv`;
- `external_reference_summary.json`;
- `external_reference_report.md`;
- the four-panel figure in `figures/`;
- per-example `comparison.json` files; and
- `provenance/sha256_manifest.txt`.

The revised Example 5 comparison is verified separately with:

```powershell
D:\anaconda3\envs\mypy\python.exe .\verify_independent_example5.py
```

This command compares the independent published-algorithm series directly
with the corrected harness MATLAB output and writes
`example5_cui/independent_run_comparison.json`.

The series-level provenance table can be regenerated with:

```powershell
D:\anaconda3\envs\mypy\python.exe .\provenance\build_external_reference_manifest.py
```

It writes `provenance/external_reference_provenance_sha256.csv`, including the
SHA-256 of each external-reference series and the corresponding harness output.

## Comparison rules

- The existing author-confirmed external outputs are treated as immutable inputs.
- The independent Example 5 output is a newly deposited numerical reference series; it does not replace or rewrite the historical files.
- The comparison script writes derived metrics and figures separately and does not overwrite the archived external series.
- The external series were not generated from the corrected harness MATLAB outputs.
- Current benchmark parameters and ground-motion file hashes were checked against the archived studies.
- Comparisons use the common sampled interval.
- No response amplitude rescaling, fitted time shift, or response sign adjustment is applied.
- Whole-mode sign alignment is permitted only for Example 1 because eigenvector sign is arbitrary.
- NRMSE is RMSE divided by the range of the external-reference series.
- The archived Example 4 acceleration series is a relative-acceleration output and is compared with the matching corrected harness MATLAB quantity.
- The archived Example 4 force series is retained, but it is not summarized because its force-output convention differs from the corresponding harness output definition.

## Distribution boundary

Downloaded program copies and book page scans are not distributed in the public validation package. Public evidence consists of source metadata, integrity hashes, archived author-confirmed outputs, the independently executed Example 5 numerical histories, and the verification scripts. The provenance records identify the original archive or project filenames; where the original generation date is not recoverable from the archived material, it is reported as unavailable rather than inferred. The series-level mapping between each external reference and its corresponding harness output is recorded in `provenance/external_reference_provenance_sha256.csv`.
