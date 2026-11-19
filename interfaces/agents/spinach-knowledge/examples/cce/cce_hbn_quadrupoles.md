# cce_hbn_quadrupoles

Syntax: `[answer,diagnostics]=cce_hbn_quadrupoles()`. The numeric answer is the projected Hahn-echo amplitude; diagnostics is the separate clexpn audit record.

The example uses standard sys/inter/bas/parameters blocks, a declarative magnetic-quantum-number transition, and `clexpn(spin_system,@hahn_echo,parameters,'labframe')`. Its field, isotopes, coordinates, coupling tensors, time grid, bath ensemble, and complete cluster graph are unchanged.

Independent dense Hamiltonian, branch propagators, and complex echo comparison live in test_cce_examples. This small model is illustrative, not a material-specific quantitative reproduction.
