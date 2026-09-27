# examples/fundamentals/pfg_test_1.m

- Signature: `pfg_test_1()`

## Purpose

Exercise the explicit gradient-pulse function grad_pulse, which uses the auxiliary-matrix formalism to compute a sample-volume integral. Background: [10.1016/j.jmr.2014.01.011](http://dx.doi.org/10.1016/j.jmr.2014.01.011).

## Physical / mathematical content

The example evolves a three-proton system under a rectangular homospoil gradient pulse and examines the trajectory by coherence order. Chemical shifts and scalar couplings are randomized for the run.

## Numerical / algorithmic content

It evaluates 100 pulse durations at steps of 2e-7 s, with gradient strength 20 G/cm, sample length 1.5 cm, and rectangular shape factor 1. The trajectory loop uses parfor.

## Implementation structure

- Set a 5.9 T, three-1H system with random scalar shifts and couplings, then construct the sphten-liouv basis and Hamiltonian.
- Build and weight the initial state by coherence-order subspaces, then call grad_pulse for each duration.
- Plot the resulting state trajectory with trajan(...,'coherence_order').
