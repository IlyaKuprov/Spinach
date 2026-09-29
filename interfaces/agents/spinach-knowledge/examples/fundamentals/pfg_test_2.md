# examples/fundamentals/pfg_test_2.m

- Signature: `pfg_test_2()`
- Source: [examples/fundamentals/pfg_test_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/pfg_test_2.m)
- Background reference: [Journal of Magnetic Resonance, DOI 10.1016/j.jmr.2014.01.011](https://doi.org/10.1016/j.jmr.2014.01.011)
- Source-stated runtime: seconds

## Purpose

Demonstrate the auxiliary-matrix algorithm for a two-gradient sandwich that selects multiple-quantum coherence. The displayed quantity is the coherence-order trajectory through the second gradient interval; the source does not calculate or report a scalar filter efficiency.

## Spin system and preparation

The no-argument function builds a three-`1H` system at 5.9 T, with randomised scalar shifts (`10*rand(1)` each) and scalar pair couplings (`20*rand(1)` each), then creates an untruncated `sphten-liouv` basis, applies the NMR assumption, and constructs the Hamiltonian. It constructs a hard `pi/2` propagator from the `1H` `Lx` operator. The initial vector is randomised, its first component is set to 1, and it is normalised. Coherence orders are obtained by summing the projection quantum numbers from `lin2lm`; the weighting vector `[0,0,0,0,0,1,0]` assigns initial weight only to the sixth enumerated coherence-order subspace, with each selected subspace normalised before weighting.

## Gradient sandwich and output

`grad_sandw` is called at 1,000 points from 0 to 2 ms. The two gradient strengths are `[20 20]` G/cm, the sample length is 2.5 cm, and both shape factors are 1 (rectangular). The durations begin as `[1 1]` ms; at each point the second duration is replaced by the current time while the first remains 1 ms. The source parallelises the calls with `parfor`, then plots `trajan(...,'coherence_order')` with a linear vertical axis.

The file supplies a plotted trajectory and a source-stated runtime estimate, but no pass/fail assertion, tolerance, or reported numerical efficiency. The random shifts and couplings are not seeded within the function, so a particular trajectory is not fixed by the source.

## Callable context

Run `pfg_test_2()` with the Spinach system, basis, Hamiltonian, operator, propagator, gradient-sandwich, coherence-analysis, and plotting functions used in the source. It accepts no arguments and produces a figure.
