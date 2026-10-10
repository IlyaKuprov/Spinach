# examples/fundamentals/quadratures/product_quadratures_1.m

- MATLAB implementation: [examples/fundamentals/quadratures/product_quadratures_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/product_quadratures_1.m)

## Purpose

Compare product-quadrature propagation schemes as the time grid is refined for an E1000B Veshtort-Griffin pulse. The source calls this an accuracy test but defines no pass/fail threshold.

## System and numerical question

The model contains 31 `1H` spins at 14.1 T. Their scalar Zeeman values are linearly spaced from -4 to 4; adjacent spins are coupled by 10 Hz. The basis uses `sphten-liouv`, the IK-2 approximation, scalar-coupling connectivity, and proximity level 1, followed by the NMR assumption. The initial state is proton `Lz`, and the control operator is proton `Lx`. The pulse duration is 0.01 s.

The comparison asks how final-state error varies with grid density for left-point and midpoint propagation, a two-point product quadrature, and a three-point product quadrature.

## Method and checks

The reference pulse has 2,000 nominal grid points; the source generates 3,999 amplitude samples and applies them in three-sample `step` calls at stride two. The benchmark grid sizes are 100, 200, ..., 1,000 points. Left-point and midpoint methods use one amplitude sample per step; the midpoint uses the average of the adjacent endpoint amplitudes. The two-point and three-point methods pass two or three adjacent Hamiltonians to `step`, respectively.

For each method the plotted error is `norm(rho_ref - rho_method) / norm(rho_ref)`. The output figure shows the pulse amplitude and error versus grid-point count; error is plotted on a logarithmic scale. No numerical tolerance or pass/fail decision is implemented.

## Output and limitations

The source provides a comparative error plot, not a tabulated result or assertion that a particular method meets a tolerance. The reference is a finer calculation using the three-point propagation family, not an analytic solution. The source header estimates runtime as seconds; no runtime was measured for this note.

[Source example](../../../../../../examples/fundamentals/quadratures/product_quadratures_1.m).
