# examples/singlet_states/m2s_example.m

- Signature: `m2s_example()`

## Purpose

Demonstrates magnetisation-to-singlet (M2S) conversion for a coupled pair of `13C` spins: the initial longitudinal state is converted by the M2S sequence and the final overlap with the singlet is displayed.

## Spin system and pulse operators

The model has two `13C` spins at `9.4 T`, scalar Zeeman values `0.03` and `-0.03`, and scalar coupling `55`. The source does not label units for those scalar values. It uses an unapproximated `sphten-liouv` basis, an NMR Hamiltonian, and carbon `Lx` and `Ly` operators as the sequence's pulse operators.

## Preparation and observable

The input is `state(spin_system,'Lz','all')`; the detector is the singlet on spins 1 and 2. The code calls `m2s` with arguments `55` and `6.0`, then displays the detector overlap with the returned state as the singlet population. These are sequence arguments as supplied by the example; this file does not define a pulse waveform, phase/amplitude/time discretisation, gradient, relaxation, or storage model. The source labels the calculation time as seconds but does not report a numerical overlap.

## Source

[examples/singlet_states/m2s_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/m2s_example.m)
