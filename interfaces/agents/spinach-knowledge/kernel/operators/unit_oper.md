# kernel/operators/unit_oper.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/unit_oper.m
Wiki: https://spindynamics.org/wiki/index.php?title=unit_oper.m

## Purpose and output

`unit_oper(spin_system)` returns a sparse identity matrix for the selected formalism. It is diagonal in, and preserves, the current basis ordering; it does not construct a generator or time propagator.

## Dimensions

Let `d=prod(spin_system.comp.mults)`. The function uses these dimensions:

- `sphten-liouv`: `size(spin_system.bas.basis,1)`.
- `zeeman-hilb` and `zeeman-wavef`: `d`.
- `zeeman-liouv`: `prod(spin_system.comp.mults.^2)`, equal to `d^2`.

For each supported formalism the returned matrix is `speye` of the specified dimension, so it acts as the identity in that representation. Any other `spin_system.bas.formalism` raises an error.
