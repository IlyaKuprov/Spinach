# cce_capabilities

Syntax: `contract=cce_capabilities(pulse_sequence)`

Explicit handle registry for free_evolution, hahn_echo, cpmg_echo, cce_eseem and bath_pumped_hahn. Declares formalism, projected support, closed-system restriction, multiplicative reduction and evaluator. Unknown handles are rejected before generator construction. The evaluator calls pulse_sequence(spin_system,parameters,H,R,K) directly.
