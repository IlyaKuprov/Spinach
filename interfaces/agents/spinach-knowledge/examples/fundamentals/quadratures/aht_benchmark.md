# examples/fundamentals/quadratures/aht_benchmark.m

- Signature: `aht_benchmark()`

## Purpose

Benchmark the accuracy of period propagators built with left-edge, midpoint, second-order Lie, and fourth-order Lie quadratures for a sinusoidally modulated RF Hamiltonian. The source describes the calculation time as seconds and notes that the second- and fourth-order methods show the same apparent accuracy for this modulation.

## Physical / mathematical content

The example sets up a single 14N spin at 14.1 T with a quadrupolar coupling and an RF field, using a rotating-frame period propagator as the reference. The source parameters include quadrupole values 1.18e6 and 0.53, RF power 2*pi*55e3/sin(theta), and RF frequency 48e3.

## Numerical / algorithmic content

A 64-slice fourth-order Lie calculation supplies the reference. For grids of 2:32 slices per RF period, each method is compared with it using norm(P_ref-P,1)/norm(P_ref,1); the error is plotted against slice count on logarithmic axes.

## Implementation structure

- Build the 14N spin system, basis, relaxation settings, and rotating-frame parameters.
- Form the reference propagator from 64 fourth-order Lie slices.
- Construct left-edge and midpoint products, then second- and fourth-order Lie products using isergen; plot their normalized errors.
