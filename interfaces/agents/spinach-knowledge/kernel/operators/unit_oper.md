# kernel/operators/unit_oper.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/unit_oper.m
Wiki: https://spindynamics.org/wiki/index.php?title=unit_oper.m

## Purpose and output

`unit_oper(spin_system)` returns a sparse identity matrix for the selected formalism. It is diagonal in, and preserves, the current basis ordering; it does not construct a generator or time propagator.

## Dimensions

The dimension is `spin_system.bas.offsets(end)`, the sum of the compiled substance dimensions. The result is `speye` of that dimension in every formalism, including a one-dimensional block for each spin-free substance.
