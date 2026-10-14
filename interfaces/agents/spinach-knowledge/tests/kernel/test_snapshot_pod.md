# Snapshot POD regression references

The regression uses independently known weighted singular values and projection residuals, equivalent covariance factors, seeded reproducibility, and caller RNG preservation. It compares full-rank and invariant-subspace complex dissipative traces against a closed-form triangular solution, and a damped proton coherence against ordinary Spinach evolution with explicit machine-precision pruning thresholds.

The test separates algebraic projection correctness from truncation accuracy, checks complex detection conjugation and Spinach's propagation sign, and rejects unsupported ranks and input shapes. Compact scaled numerical references target 1e-12 maximum absolute error; they do not certify arbitrary large or ill-conditioned reductions.

At scale 1e-200, exact SVD and randomised sketches with zero/two power iterations check scalar singular values and rank-two truncated subspaces using scaled or relative errors. Rounded-zero energy is accepted while zero data, unsupported rank, and energy overflow remain errors; caller RNG is checked on acceptance and rejection.

Mixed sparse/dense weighted blocks with a zero-weight block exercise genuinely truncated sketches with zero and two power iterations against an explicit concatenated reference. Both exact and randomised methods compare sparse snapshots with dense equivalents and reject stored NaN/Inf and all-zero data without changing caller RNG. Sparse projection accepts a zero Liouvillian term and rejects stored NaN/Inf entries.

A zero-step trajectory is checked against the independent initial complex overlap. Zero-duration traces are independent of generator assembly, including a case whose affine assembly would overflow. Zero-weight column exclusion is checked with active data narrower than the requested sketch, padded singular values, rank rejection, and unchanged caller RNG. Invalid zero-weight blocks are still rejected.

Real truncated sketches with zero and two power iterations match an independent real Gaussian reference, retain real basis storage despite a zero-weight complex block, and preserve RNG and seed repeatability. Dense NaN/Inf terms are rejected alongside the sparse validation cases.

Sparse coordinate-selection preparation checks reduced storage, the identity projection, complex Hermitian detection, and a nonzero-duration analytic phase trajectory. Sparse reduced models cover zero-step and zero-duration queries, a zero term, and rejection of stored NaN/Inf before initial-only returns. Exact all-zero coefficients are tested for dense and large sparse models with complex initial overlaps; stored NaN/Inf terms remain invalid for that shortcut.

Constant-signal cases cover full and sparse zero coefficient columns, zero sparse terms, exact affine cancellation, and zero states/detectors. A two-coordinate mixing reference ensures a zero initial overlap can evolve into a nonzero signal. Stored NaN/Inf remain invalid with zero coefficients.

Large oversampling retains the prescribed leading subspace and all compressed singular values. Exact sparse right/left nullspace fixtures include non-normal dynamics and a nonstationary state observed by a stationary detector. A tiny nonzero phase is not discarded as stationary.

Underflowing state and detector actions retain their finite independent phase reference. A cancellation-based nullspace remains supported without the structural shortcut.

Dense/sparse real and complex state projection overflow is rejected during preparation; large finite projected initial and detector coordinates remain supported.

Dense/sparse real and complex initial-overlap overflow is rejected before zero-step, zero-duration and zero-coefficient returns, and ordinary positive-time requests. Non-cancelling product overflow is also rejected. Large cancelling overlaps with finite intermediate products remain supported.

Time-scaled affine assembly is checked with a detector-visible tiny mixing phase, coefficient/time underflow and overflow, an overflowing unscaled generator with a finite damped exponent, and explicit rejection of lost nonzero scaled entries. Later state and real/complex detector overflow are rejected, with independently normalised finite amplified controls. Dense and sparse terms exercise the same regimes.

Dense/sparse real and imaginary fixtures check subnormal product-order fallback, affine cancellation before time scaling with a nonzero residual, and common state/detector structural stationarity despite unreachable overflow. Mixed-side stationarity and lost tiny entries remain adversarial controls.
