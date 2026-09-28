# kernel/ppower.m

- Signature: `P=ppower(spin_system,P,N)`

## Purpose

Raises the propagator matrix `P` to a non-negative integer power `N`. The implementation processes the binary digits of `N`, multiplying only the powers corresponding to set bits.

## Mathematical and numerical behaviour

For `N=0`, the result is the identity matrix; for `N=1`, the input propagator is returned unchanged. For larger powers, repeated squaring builds the required powers of `P`, and an accumulator combines the terms selected by the binary expansion of `N`. Each accumulator multiplication and required squaring is passed through `clean_up` with `spin_system.tols.prop_chop`. Sparse inputs use a sparse identity; otherwise the identity is created with the input matrix's numeric type.

## Parameters / inputs

- `spin_system` — Spinach spin system structure; it must contain a non-negative real scalar `spin_system.tols.prop_chop`.
- `P` — square numeric propagator matrix.
- `N` — non-negative real integer power, supplied as a scalar. Non-integer numeric inputs are limited to `flintmax` before conversion to `uint64`.

## Output

- `P` — propagator matrix raised to the power `N`; the zero-power result is an identity matrix.

The routine checks the spin-system tolerance, matrix shape, and power argument before computing (except that the implementation's `N=1` shortcut follows those checks).

## Reference

[Spin Dynamics Wiki: `ppower.m`](https://spindynamics.org/wiki/index.php?title=ppower.m)
