# kernel/average.m

- Signature: `H=average(spin_system,Hp,H0,Hm,omega,theory)`

## Purpose

Builds an average Hamiltonian from the positive-, zero-, and negative-frequency components of a Hamiltonian in a Zeeman-interaction rotating frame.

## Physical / mathematical content

- `Hp`, `H0`, and `Hm` are the rotating-frame Hamiltonian components with frequencies `+omega`, `0`, and `-omega`, respectively. The selected theory determines the terms included in the average Hamiltonian.
- The Waugh average-Hamiltonian (AH) options retain terms through first, second, or third order. The Krylov-Bogolyubov (KB) options provide the corresponding first-, second-, or third-order theory for DNP experiments.
- `matrix_log` is described by the source as an exact, expensive dense-matrix method. Its implementation forms a one-period product integral using 16 time slices and fourth-order Lie quadrature, then applies a matrix logarithm.

## Parameters / inputs

- `spin_system` - spin-system structure used for reporting diagnostics.
- `Hp`, `H0`, `Hm` - numeric matrices of equal size, as defined above.
- `omega` - non-zero finite real rotating-frame frequency in rad/s.
- `theory` - character string selecting one of:
  - `ah_first_order`, `ah_second_order`, `ah_third_order` - first-, second-, or third-order Waugh theory.
  - `kb_first_order`, `kb_second_order`, `kb_third_order` - first-, second-, or third-order Krylov-Bogolyubov theory for DNP experiments.
  - `matrix_log` - the product-integral/matrix-log method described above.

## Outputs

- `H` - average Hamiltonian.

## References

Krylov-Bogolyubov averaging for DNP systems is described in:

- [10.1039/C2CP23233B](http://dx.doi.org/10.1039/C2CP23233B)
- [10.1007/s00723-012-0367-0](http://dx.doi.org/10.1007/s00723-012-0367-0)

[Spin Dynamics Wiki: average.m](https://spindynamics.org/wiki/index.php?title=average.m)
