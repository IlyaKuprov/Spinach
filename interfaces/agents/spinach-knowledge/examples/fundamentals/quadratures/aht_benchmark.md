# examples/fundamentals/quadratures/aht_benchmark.m

- Signature: `aht_benchmark()`
- Source: [examples/fundamentals/quadratures/aht_benchmark.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/aht_benchmark.m)
- Source-stated runtime: seconds

## Purpose

Benchmark how several slice-based propagator constructions approximate one period of a periodically modulated Hamiltonian. It compares left-edge and midpoint sampling with second- and fourth-order Lie-group quadratures against a finer fourth-order reference. The source notes that the modulation is purely sinusoidal, so second- and fourth-order integrators show the same apparent accuracy in this example; this is the source's stated expectation, not a pass criterion.

## Model and rotating-frame setup

The no-argument example uses one `14N` spin at 14.1 T, an `eeqq2nqi` quadrupolar coupling built from `1.18e6`, asymmetry `0.53`, spin `1`, and Euler angles `[0 0 0]`, plus scalar Zeeman value `32.4`. The basis is `sphten-liouv` without truncation. Damping relaxation is configured with diagonal retention, zero equilibrium, and damping rate 300; the `krylov` and `trajlevel` options are disabled. The pulse setup uses the magic angle `theta=atan(sqrt(2))`, single-crystal orientation `[sqrt(2/3) 0 sqrt(1/3)]`, rank 6, rate `-19840`, RF power `2*pi*55e3/sin(theta)`, and RF frequency `48e3`. The Hamiltonian and parameters are extracted through `singlerot(...,@impound,...,'qnmr')`; the source then forms the overtone-offset frequency and two complex-conjugate sinusoidal RF Hamiltonian components.

## Numerical comparison and output

The reference product uses 64 slices and `isergen(HL,HM,HR,slice_dur)` at left-edge, midpoint, and right-edge Hamiltonians. For each slice count from 2 through 32, the script forms four products: left-edge constant Hamiltonians, midpoint constant Hamiltonians, second-order Lie quadrature from the left and right values, and fourth-order Lie quadrature using left, midpoint, and right values. It plots each relative propagator error, `norm(P_ref-P,1)/norm(P_ref,1)`, against slice count on logarithmic axes.

The source defines a numerical comparison, not a pass/fail test: it sets no acceptance threshold and contains no assertion or printed table of measured error values. The reference is itself a finite-slice fourth-order calculation, not an analytic exact propagator.

## Callable context

Run `aht_benchmark()` with the Spinach single-rotation, propagator, generator, and plotting functions called by the source. It takes no arguments and produces the error plot.
