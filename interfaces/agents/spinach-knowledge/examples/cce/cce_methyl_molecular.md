# cce_methyl_molecular

Syntax: `[answer,diagnostics]=cce_methyl_molecular()`. The numeric answer is the projected Hahn-echo amplitude; diagnostics is the separate clexpn audit record.

The example uses standard sys/inter/bas/parameters blocks, a declarative magnetic-quantum-number transition, and `clexpn(spin_system,@hahn_echo,parameters,'labframe')`. Its field, isotopes, coordinates, coupling tensors, time grid, bath ensemble, and complete cluster graph are unchanged.

Independent dense Hamiltonian, branch propagators, and complex echo comparison live in test_cce_examples. This small model is illustrative, not a material-specific quantitative reproduction.

The disjoint expansion uses one three-proton group. The exchange terms are coherent pair-exchange interactions; they are not a stochastic methyl-jump model.
