# Snapshot POD regression references

The regression uses independently known weighted singular values and projection residuals, equivalent covariance factors, seeded reproducibility, and caller RNG preservation. It compares full-rank and invariant-subspace complex dissipative traces against a closed-form triangular solution, and a damped proton coherence against ordinary Spinach evolution with explicit machine-precision pruning thresholds.

The test separates algebraic projection correctness from truncation accuracy, checks complex detection conjugation and Spinach's propagation sign, and rejects unsupported ranks and input shapes. Compact scaled numerical references target 1e-12 maximum absolute error; they do not certify arbitrary large or ill-conditioned reductions.

At scale 1e-200, exact SVD and randomised sketches with zero/two power iterations check scalar singular values and rank-two truncated subspaces using scaled or relative errors. Rounded-zero energy is accepted while zero data, unsupported rank, remain errors; caller RNG is checked on acceptance and rejection.

Mixed sparse/dense weighted blocks with a zero-weight block exercise genuinely truncated sketches with zero and two power iterations against an explicit concatenated reference. Both exact and randomised methods compare sparse snapshots with dense equivalents and reject stored NaN/Inf and all-zero data without changing caller RNG. Sparse projection accepts a zero Liouvillian term and rejects stored NaN/Inf entries.

A zero-step trajectory is checked against the independent initial complex overlap. Zero-duration traces are independent of generator assembly, including a case whose affine assembly would overflow. Zero-weight column exclusion is checked with active data narrower than the requested sketch, padded singular values, rank rejection, and unchanged caller RNG. Invalid zero-weight blocks are still rejected.

Real truncated sketches with zero and two power iterations match an independent real Gaussian reference, retain real basis storage despite a zero-weight complex block, and preserve RNG and seed repeatability. Dense NaN/Inf terms are rejected alongside the sparse validation cases.

Sparse coordinate-selection preparation checks reduced storage, the identity projection, complex Hermitian detection, and a nonzero-duration analytic phase trajectory. Sparse reduced models cover zero-step and zero-duration queries, a zero term, and rejection of stored NaN/Inf before initial-only returns. Exact all-zero coefficients are tested for dense and large sparse models with complex initial overlaps; stored NaN/Inf terms remain invalid for that shortcut.

Constant-signal cases cover full and sparse zero coefficient columns, zero sparse terms, exact affine cancellation, and zero states/detectors. A two-coordinate mixing reference ensures a zero initial overlap can evolve into a nonzero signal. Stored NaN/Inf remain invalid with zero coefficients.

Large oversampling retains the prescribed leading subspace and all compressed singular values. Exact sparse right/left nullspace fixtures include non-normal dynamics and a nonstationary state observed by a stationary detector. A tiny nonzero phase is not discarded as stationary.

Underflowing state and detector actions retain their finite independent phase reference. A cancellation-based nullspace remains supported without the structural shortcut.

Dense/sparse real and complex state projection overflow is rejected during preparation; large finite projected initial and detector coordinates remain supported.

Dense/sparse real and complex initial overlaps use bounded detection before physical rescaling: opposing large products can cancel to finite signals in constant and ordinary queries. Non-cancelling unrepresentable physical signals remain errors.

Time-scaled affine assembly is checked with a detector-visible tiny mixing phase, coefficient/time underflow and overflow, an overflowing unscaled generator with a finite damped exponent, and explicit rejection of lost nonzero scaled entries. Finite detector-scaled amplification beyond physical state range is supported; unrepresentable real/complex detector signals remain errors, with independently normalised finite controls. Dense and sparse terms exercise the same regimes.

Exponent-range boundaries include subnormal coefficient/time products, both finite factor orderings, and exact affine cancellation before an otherwise overflowing time scale. Unreachable state/detector entries do not require representable scaling. Two terms stationary on different sides compose a detector-visible nilpotent chain, preventing a per-term mixed-side shortcut.

Dense/sparse real and imaginary fixtures check subnormal product-order fallback, affine cancellation before time scaling with a nonzero residual, and common state/detector structural stationarity despite unreachable overflow. Mixed-side stationarity and lost tiny entries remain adversarial controls.

Power-of-two-normalised exact and randomised training checks scale covariance at1e-200,1 and1e200, with independent singular values and residual ratios. Physical squared energy is a diagnostic, not a rank decision. Caller RNG remains preserved.

Compensating block amplitudes and weights have independently known singular values and reconstruction residuals. Physical diagnostic overflow retains a finite scaled frame; unsupported component loss is rejected without RNG changes. Native complex scalar SVD bases are accepted by the dimension-derived Gram budget. A positively amplified non-normal triangular family matches an independently scaled analytic semigroup beyond physical state range.

Streamed dense-term projection is checked against independently selected coordinates of a1000-by1000 complex operator; projected state and detector selections are exact. Sparse prepared storage remains a maintained acceptance condition. Non-double time and integer-count inputs are rejected explicitly.

Additional references cover scaled-product rejection, exponent-aware real and complex overlap, disconnected unstable modes, and directed source-to-detector paths against independent triangular formulas.

Additional real/imaginary tiny-basis fixtures reject product loss in both Galerkin stages. Constant-query detector fixtures now require the finite exponent-aware overlap instead of rejecting it. Exact complex-component, projection and detection cancellations remain independent controls.

Dense real/complex invariant-mode references cover short and composed small-step actions despite an unobserved unstable mode. Adjoint references verify the propagation sign and complex Hermitian detection. Genuine and eps-sized observed growing modes are adversarial controls; neither is truncated as invariant. Native action RNG is checked on acceptance and error.

Complete-support dense/sparse normal dynamics use an independent rank-one exponential formula, with real/complex detection and one/two affine terms. These references validate graph bypass and callback range caching without inferring accuracy from their implementation.

Dense-basis diagonal projection uses independent parity-sum formulas for real/complex coordinates with sparse/dense terms. Incomplete dense supports and disconnected controls compare implicit directed traversal against independent full exponentials, with real/complex detection and both storage types.

Real-storage controls compare dense/sparse coordinate projections with exact selected entries, and retain genuinely complex terms and states.

Independent rank-one relaxation exponentials and two-mode triangular formulas exercise separated diagonals, real and complex coefficients, Hermitian detection, and exchanged source/detector roles in dense and sparse storage.

Scalar complex-phase controls preserve the selected adjoint conjugation in the closed-form action, including exchanged state and detector roles.

Connected tridiagonal relaxation fixtures compare dense and rank-8192 sparse models with an independent alternating-eigenmode exponential, complex coefficients, and exchanged source/detector roles.

Affine connecting-path cancellations exclude an unobserved growing mode before restoring its physical scale, including a sum whose irrelevant diagonal would overflow. Dense/sparse and real/imaginary detectors use an independent scalar phase; residual connecting paths remain explicit range errors.

All six orders of a three-term affine family preserve a weak nilpotent coupling between cancelling large entries, with dense/sparse terms and real/imaginary source phases. A five-term deeper cancellation must explicitly reject an uncertified zero; an exactly accumulated reordering recovers its independent nilpotent signal.
