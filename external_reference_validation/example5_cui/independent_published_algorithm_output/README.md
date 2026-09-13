# Example 5 independent published-algorithm output

These response histories were obtained by independently transcribing and
executing the MATLAB algorithm published in Chapter 11 of Cui, Shen, and Yang
(2022), together with the required Appendix 1 and Appendix 6 functions, in
MATLAB R2025b. The ground-motion record was normalized to PGA = 0.40 g only to
match the benchmark input convention used by the corrected harness MATLAB
output.

The files contain numerical response histories in SI units:

- `top_displacement_m.txt`: time (s), top displacement (m);
- `top_relative_acceleration_mps2.txt`: time (s), top relative acceleration (m/s^2);
- `isolator_hysteresis_m_N.txt`: isolation displacement (m), restoring force (N).

The source listing and page scans are copyrighted and are not redistributed in
this repository. The output files are deposited so that the revised Example 5
comparison can be independently audited. No response-history amplitude
rescaling, fitted time shift, or response sign adjustment was applied after the
published algorithm was run.
