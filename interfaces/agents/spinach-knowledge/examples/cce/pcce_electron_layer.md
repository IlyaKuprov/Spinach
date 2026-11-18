# examples/cce/pcce_electron_layer.m

Syntax: `[answer,diagnostics]=pcce_electron_layer()`. The numeric answer is the projected Hahn-echo amplitude; diagnostics is the separate clexpn audit record.

The example uses standard sys/inter/bas/parameters blocks, a declarative magnetic-quantum-number transition, and `clexpn(spin_system,@hahn_echo,parameters,'labframe')`. The observed electron is coupled to three electron bath spins with a mixed-product bath state and complete graph. Partitioned order two retains the strongly coupled pair `[2 3]` as one group and spin 4 as the other. Explicit coupling tensors, rather than coordinate-derived couplings, define this illustrative planar model; the 41 total echo times run from zero to 200 microseconds.

The example plots the real and imaginary parts of the normalised complex echo. `tests/kernel/test_cce_examples.m` compares it with a separate unpartitioned full-order expansion and an independently assembled dense projected-bath calculation, and retains the reference physical arrays. This finite model is not a quantitative material reproduction.

[MATLAB source](../../../../../examples/cce/pcce_electron_layer.m) · [Independent reference tests](../../../../../tests/kernel/test_cce_examples.m)
