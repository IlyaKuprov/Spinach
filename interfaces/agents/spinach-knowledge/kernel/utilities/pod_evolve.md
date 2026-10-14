# Reduced observable trajectories for an explicit POD model

`pod_evolve` evaluates a prepared projected affine model with an explicit coefficient vector and time grid. The initial complex observable is included. The evaluator assembles the time-scaled affine generator before MATLAB's fixed dense exponential, using Spinach's exp(-1i*L*t) sign, with L in rad/s and time in seconds; no physical spin_system or propagation screening policy is modified.

This evaluator does not train, retune rank, subtract a reference, or process a spectrum. Its intended use is an explicitly validated fixed family for which projection costs can be amortised. Generator assembly and reduction errors remain the caller's responsibility; no observable error bound follows from orthogonality.

Zero-step, zero-duration, and exactly all-zero coefficient requests return the initial observable, repeated as needed, after validation, without assembling a generator or computing a matrix exponential. The coefficient shortcut uses exact zeros; it does not infer cancellation between nonzero terms.

Sparse model terms and state vectors use finiteness checks that preserve implicit zeros. This applies before initial-only returns: a sparse zero term is accepted, while stored NaN/Inf entries are rejected. Assembly omits zero coefficients and retains a sparse generator when active terms are sparse; nonzero propagation expands the generator for the dense exponential.

All-zero coefficients, an exactly zero initial state or detector, and an exactly zero assembled generator return constant signals without an exponential. These are exact algebraic cases, not tolerance screening; zero initial overlap alone does not imply a constant trajectory. Invalid model terms are rejected even when coefficients are zero.

Zero generator columns on the initial-state support, or zero rows on the detector support, certify a constant signal before dense expansion. This structural check uses stored zeros, not computed action magnitudes: a floating-point action can underflow despite a finite observable phase. Nullspaces that require cancellation remain supported by the general exponential path.

The computed initial overlap must be finite before any initial-only or constant-signal return. Finite coordinates do not guarantee finite intermediate products: large opposing products can overflow before cancelling. Such overlaps are rejected explicitly; no overflow-safe rescaling is inferred.

Coefficient/time products are used when normal and finite; otherwise each term is multiplied by its coefficient before the time step, avoiding a lost intermediate scalar scale. Nonzero term entries lost to time-scaling underflow are rejected rather than treated as exact structural zeros. Computed propagated states and observables must remain finite; no scaled-state propagation is inferred.

Common zero state-support columns or detector-support rows across every active term certify stationarity before scaling unreachable entries; conditions from opposite sides cannot be mixed across terms. Nonfinite term-first products retry a finite nonzero coefficient/time scale. Nonfinite separately scaled assembly retries the unscaled affine sum before time scaling, with finite products and nonzero-entry representability still required. These bounded fallbacks do not guarantee support for every mathematically finite exponent.
