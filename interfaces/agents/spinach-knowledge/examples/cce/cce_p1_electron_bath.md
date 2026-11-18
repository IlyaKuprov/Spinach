# cce_p1_electron_bath

Syntax: `[answer,diagnostics]=cce_p1_electron_bath()`. The numeric answer is the projected Hahn-echo amplitude; diagnostics is the separate clexpn audit record.

The example uses standard sys/inter/bas/parameters blocks, a declarative magnetic-quantum-number transition, and `clexpn(spin_system,@hahn_echo,parameters,'labframe')`. Physical tensors, elapsed-time grid, graph, retained order, and bath state are unchanged.

Independent analytic, source-pinned, or dense checks and frozen physical arrays are retained in `tests/kernel/test_cce_examples.m`.
