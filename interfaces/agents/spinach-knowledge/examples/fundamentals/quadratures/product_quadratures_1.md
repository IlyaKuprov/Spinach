# examples/fundamentals/quadratures/product_quadratures_1.m

- Signature: `product_quadratures_1()`

## Purpose

Tests the accuracy of Lie-group product quadratures as the time-grid spacing changes during an E1000B Veshtort–Griffin pulse.

## Physical / mathematical content

- The model contains 31 `1H` spins at 14.1 T, with Zeeman shifts linearly spaced from -4 to 4 and 10 Hz scalar couplings between adjacent spins. The basis is `sphten-liouv` with the IK-2 approximation, scalar-coupling connectivity, and proximity level 1.
- The pulse lasts 10 ms. The initial density operator is `Lz` on `1H`, and the control operator is `Lx` on `1H`.

## Numerical / algorithmic content

- A 2000-point, three-Hamiltonian-step calculation supplies the reference. The benchmark compares left-point, midpoint, two-point Lie-group, and three-point Lie-group propagation (LP, MP, LG-2, LG-4) at 100, 200, ..., 1000 grid points, recording relative state errors.

## Implementation structure

- Generates the E1000B pulse with `vg_pulse`, forms the reference, and runs the benchmark cases in a `parfor` loop. The figure shows pulse amplitude and relative error versus time-grid size on log-log axes.
- The source estimates a runtime of seconds.
