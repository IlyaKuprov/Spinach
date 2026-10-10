# examples/fundamentals/pfg_test_1.m

- Signature: `pfg_test_1()`
- Source: [examples/fundamentals/pfg_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/pfg_test_1.m)
- Background reference: [Journal of Magnetic Resonance, DOI 10.1016/j.jmr.2014.01.011](https://doi.org/10.1016/j.jmr.2014.01.011)

## Purpose

Demonstrate an explicit gradient pulse calculation using the auxiliary-matrix formalism to obtain the sample-volume-integrated state trajectory. The observable displayed is how coherence-order content evolves during a homospoil gradient; the source does not reduce it to a reported scalar filter efficiency.

## Spin system and initial state

The no-argument example creates three `1H` spins at a 5.9 T magnet, with independently randomised scalar shifts (`10*rand(1)` per spin) and scalar pair couplings (`20*rand(1)` for each pair). It uses the `sphten-liouv` basis with `approximation='none'` and applies the NMR assumption. A normalised random state vector is used as the initial state, with its first element set to 1 before normalisation. Basis projection quantum numbers from `lin2lm` are summed to assign each basis state a coherence order. The state is normalised separately within each represented coherence-order subspace and those subspaces are weighted linearly from 0.1 to 0.9.

## Pulse calculation and output

The script calls `grad_pulse(spin_system,L,rho,gradient_strength,sample_length,duration,gradient_shape_factor)` for 100 durations from zero to 19.8 microseconds in increments of 0.2 microseconds. It sets the gradient to 20 G/cm, the sample length to 1.5 cm, and the shape factor to 1 (rectangular). The calls are parallelised with `parfor`. `trajan(...,'coherence_order')` displays the resulting trajectory with a linear vertical axis.

These are source inputs and plotted output, not a reported numerical validation result: the function has no assertion, tolerance, or scalar result. Because shifts and couplings are randomised without setting the RNG state, a run need not reproduce a particular trajectory.

## Callable context

Run `pfg_test_1()` with the Spinach system, basis, Hamiltonian, gradient-pulse, coherence-analysis, and plotting functions used in the source; it accepts no arguments and produces a figure. No external experimental data are read.

Coherence labels come from the single-substance descriptor `bas.basis{1}`, not the descriptor-cell container.
