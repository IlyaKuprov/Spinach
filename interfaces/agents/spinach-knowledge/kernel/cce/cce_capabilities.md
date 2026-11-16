# cce_capabilities

Syntax: `contract=cce_capabilities(pulse_sequence)`

Explicit handle registry for free_evolution, hahn_echo, cpmg_echo, cce_eseem and bath_pumped_hahn. Declares formalism, projected support, closed-system restriction, multiplicative reduction and evaluator. Unknown handles are rejected before generator construction. The evaluator calls pulse_sequence(spin_system,parameters,H,R,K) directly.

## Native dissipative, noise and ensemble contract

Explicit function-handle registry, version 2. Free evolution and Hahn echo admit dissipation; free evolution also admits elementwise density output. Noise correlation is additive, requests both outputs and retains moment/first-moment metadata. Unknown handles are rejected. The evaluator binds the physical configuration and invokes the ordinary `(spin_system,parameters,H,R,K)` contract.

Hahn echo also supports generalized `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.
