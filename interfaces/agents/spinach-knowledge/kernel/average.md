# kernel/average.m

Canonical implementation: [kernel/average.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/average.m)

- Signature: `H=average(spin_system,Hp,H0,Hm,omega,theory)`

## Contract

The three Hamiltonian components are numeric matrices, each must be a matrix and all three must have the same dimensions. The routine does not check their dimensions against `spin_system`; its validator also does not explicitly require square matrices. `omega` must be a finite, non-zero real scalar, in radians per second. The diagnostics display `omega/(2*pi)` in Hz. The source does not state the unit convention for `Hp`, `H0`, or `Hm`. `theory` is a character string naming one of the choices below. The result is assembled from these same-size operators; `spin_system` is used for reports.

## Theory choices

- `ah_first_order`, `ah_second_order`, and `ah_third_order`: Waugh average-Hamiltonian expressions through the named order.
- `kb_first_order`, `kb_second_order`, and `kb_third_order`: Krylov-Bogolyubov expressions through the named order; the source documents these options for DNP experiments. For example, first order is `H=H0+(Hp*Hm-Hm*Hp)/omega`. At second order the source adds `H2=-(Hp*(Hm*H0-H0*Hm)+Hm*(Hp*H0-H0*Hp))/omega^2`; the third-order implementation adds its explicit third-order expression.
- `matrix_log`: the source calls this an exact but expensive dense-matrix algorithm. It builds one rotating-frame period `2*pi/omega` using 16 intervals and fourth-order Lie quadrature, forms the product integral, and applies `logm`. This is not one of the truncated AH/KB options.

Unknown theory strings raise an error. The routine reports input/operator 1-norms, the frequency in Hz, and output dimension, nonzero count, density, 1-norm, and sparsity.

## Parameters / output

- `spin_system` - Spinach structure passed to diagnostic reporting.
- `Hp`, `H0`, `Hm` - equal-size numeric matrix components at frequencies `+omega`, zero, and `-omega` in the rotating frame.
- `omega` - rotating-frame angular frequency in rad/s; it must be finite, real, scalar, and non-zero.
- `theory` - one of the six named AH/KB orders or `matrix_log`.
- `H` - the assembled average Hamiltonian.

## References

Krylov-Bogolyubov averaging for DNP systems is described in:

- [10.1039/C2CP23233B](http://dx.doi.org/10.1039/C2CP23233B)
- [10.1007/s00723-012-0367-0](http://dx.doi.org/10.1007/s00723-012-0367-0)

[Spin Dynamics Wiki: average.m](https://spindynamics.org/wiki/index.php?title=average.m)
