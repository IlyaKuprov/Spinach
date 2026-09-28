# examples/fundamentals/pfg_test_2.m

- Signature: `pfg_test_2()`

## Purpose

Demonstrate the auxiliary-matrix algorithm for generating a gradient-sandwich multiple-quantum filter. Background: [10.1016/j.jmr.2014.01.011](http://dx.doi.org/10.1016/j.jmr.2014.01.011). The source labels the calculation time as seconds.

## Physical / mathematical content

A three-proton system is prepared with a pi/2 1H pulse. The initial state is weighted to select a coherence-order subspace, then evolved through a two-pulse gradient sandwich; the output trajectory is analysed by coherence order.

## Numerical / algorithmic content

The source uses two rectangular gradients of 20 G/cm, a 2.5 cm sample, and two 1 ms durations. It samples 1000 points over the total 2 ms trajectory and evaluates them with parfor.

## Implementation structure

- Specify random scalar shifts and couplings, create the sphten-liouv basis, and build the Hamiltonian.
- Construct a hard pi/2 pulse propagator and weight the initial state by coherence order.
- Call grad_sandw across the time axis and display the coherence-order trajectory with trajan.
