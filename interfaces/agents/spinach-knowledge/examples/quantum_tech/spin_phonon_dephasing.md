# examples/quantum_tech/spin_phonon_dephasing.m

- Signature: `spin_phonon_dephasing()`

## Purpose

Longitudinal spin–phonon coupling modulates spin coherence and displaces a quantised vibrational mode conditionally on the spin. This is a minimal Weyl-algebra version of strain-modulated spin Hamiltonians used for NV-centre and molecular spin–phonon dynamics. Both observables come from one trajectory: because the drift commutes with spin projection, the coherence and population sectors of the initial condition evolve and are detected without mixing. Calculation time: seconds.

## Physical / mathematical content

- The model couples an electron spin longitudinally to a seven-level phonon mode. It follows the spin coherence and the oscillator quadrature `a+a^+`, exposing their coupled dynamics without mixing the relevant sectors.

## Numerical / algorithmic content

- In the `zeeman-hilb` formalism with no basis approximation, Spinach propagates a mixed initial state through the `spin-phonon` device context. The sequence has 501 points and checks that both coherence modulation and oscillator displacement are visible.

## Implementation structure

- The phonon frequency is `20e6`, longitudinal coupling is `4e6*sqrt(2)`, and sequence sweep is `5e8`. The initial state is the sum of spin `Lx` and `ZL2` terms, each with phonon state `BL1`; observables are spin `Lx` and the phonon quadrature `C+A` for mode 2.
