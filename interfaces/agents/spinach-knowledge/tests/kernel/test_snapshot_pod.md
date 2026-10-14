# Snapshot POD regression references

The regression uses independently known weighted singular values and projection residuals, equivalent covariance factors, seeded reproducibility, and caller RNG preservation. It compares full-rank and invariant-subspace complex dissipative traces against a closed-form triangular solution, and a damped proton coherence against ordinary Spinach evolution with explicit machine-precision pruning thresholds.

The test separates algebraic projection correctness from truncation accuracy, checks complex detection conjugation and Spinach's propagation sign, and rejects unsupported ranks and input shapes. Compact scaled numerical references target 1e-12 maximum absolute error; they do not certify arbitrary large or ill-conditioned reductions.

Mixed sparse/dense weighted blocks with a zero-weight block exercise genuinely truncated sketches with zero and two power iterations against an explicit concatenated reference. Both exact and randomised methods compare sparse snapshots with dense equivalents and reject stored NaN/Inf and all-zero data without changing caller RNG. Sparse projection accepts a zero Liouvillian term and rejects stored NaN/Inf entries.
