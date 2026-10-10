# examples/singlet_states/s2m_example.m

- Signature: `s2m_example()`

## Purpose

Demonstrates singlet-to-magnetisation (S2M) conversion for a coupled pair of `13C` spins: the initial singlet is converted by the S2M sequence and the final longitudinal-state overlap is displayed.

## Spin system and pulse operators

The model has two `13C` spins at `9.4 T`, scalar Zeeman values `0.03` and `-0.03`, and scalar coupling `55`. The source does not label units for those scalar values. It uses an unapproximated `sphten-liouv` basis, an NMR Hamiltonian, and carbon `Lx` and `Ly` operators as the sequence's pulse operators.

## Preparation and observable

The input is the singlet on spins 1 and 2; the detector is `coil_state(spin_system,'Lz','all','exact')`. The code calls `s2m` with arguments `55` and `6.0`, then displays the detector overlap with the returned state as longitudinal magnetisation. These are sequence arguments as supplied by the example; this file does not define a pulse waveform, phase/amplitude/time discretisation, gradient, relaxation, or storage model. The source labels the calculation time as seconds but does not report a numerical overlap.

## Source

[examples/singlet_states/s2m_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/s2m_example.m)
