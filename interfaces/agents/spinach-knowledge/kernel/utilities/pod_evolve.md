# Reduced observable trajectories for an explicit POD model

`pod_evolve` evaluates a prepared projected affine model with an explicit coefficient vector and time grid. The initial complex observable is included. MATLAB's fixed dense exponential uses Spinach's exp(-1i*L*t) sign, with L in rad/s and time in seconds; no physical spin_system or propagation screening policy is modified.

This evaluator does not train, retune rank, subtract a reference, or process a spectrum. Its intended use is an explicitly validated fixed family for which projection costs can be amortised. Generator assembly and reduction errors remain the caller's responsibility; no observable error bound follows from orthogonality.

Zero-step, zero-duration, and exactly all-zero coefficient requests return the initial observable, repeated as needed, after validation, without assembling a generator or computing a matrix exponential. The coefficient shortcut uses exact zeros; it does not infer cancellation between nonzero terms.

Sparse model terms and state vectors use finiteness checks that preserve implicit zeros. This applies before initial-only returns: a sparse zero term is accepted, while stored NaN/Inf entries are rejected. Nonzero-duration propagation still explicitly forms a dense reduced generator and exponential.
