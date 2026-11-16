# cce_capabilities

Syntax: `contract=cce_capabilities(pulse_sequence)`

Explicit handle registry for free_evolution, hahn_echo, cpmg_echo, cce_eseem and bath_pumped_hahn. Declares formalism, projected support, closed-system restriction, multiplicative reduction and evaluator. Unknown handles are rejected before generator construction. The evaluator calls pulse_sequence(spin_system,parameters,H,R,K) directly.

## Native dissipative, noise and ensemble contract

Explicit function-handle registry, version 2. Free evolution and Hahn echo admit dissipation; free evolution also admits elementwise density output. Noise correlation is additive, requests both outputs and retains moment/first-moment metadata. Unknown handles are rejected. The evaluator binds the physical configuration and invokes the ordinary `(spin_system,parameters,H,R,K)` contract.

Hahn echo also supports generalized `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.
The free-evolution evaluator forwards experiment-owned coherent absolute error bounds and spectral residuals to raw records. Bounds are numerical estimates, not cluster-truncation estimates; they are not replaced by zeros.

The registered Hahn-echo evaluator transports the optional spectral/error auxiliaries exactly as for coherent free evolution. Experiments without these auxiliaries retain their existing contract.

finite_drive is registered for conventional Hilbert and Liouville generators, signal/observable/density output and density reduction; projected branches are rejected. Its auxiliary numerical error estimate is transported to the reducer.

Bath-pumped Hahn evaluation uses `clexpn(spin_system,@bath_pumped_hahn,parameters,assumptions)`. The five-argument experiment optionally returns `[E,aux]`; coherent projected propagation supplies spectral residuals, Hamiltonian round-off bounds and signal error estimates to the native reducer. The pump selects canonical axial Zeeman frequencies, with explicit Hz or rad/s units, an inclusive nonnegative halfwidth, a fraction in [0,1], logical allow_empty, and optional distinct physical bath sites or isotope-label cells. No-field product-density ensembles allow selected spins outside a retained cluster; static exterior mean-field pumping remains explicitly unsupported.

## Generalized Hahn evidence

The registered evaluator retains aux.density when a unitary generalized Hahn experiment supplies it, alongside its spectral residuals, Hamiltonian assembly bound, backend and absolute response bounds. This is per-cluster physical evidence, not a legacy public result object.

The hahn_echo contract version is 3; other experiment versions are unchanged. This separates its retained density/numerical-evidence checkpoint schema from earlier generalized raw records.
