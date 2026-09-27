# examples/quantum_tech/spin_cavity_purcell_effect.m

- Signature: `spin_cavity_purcell_effect()`

## Purpose

Demonstrates cavity-induced spin relaxation in the EPR Purcell regime: coherent Jaynes–Cummings exchange combined with rapid cavity damping in Liouville space relaxes the spin excitation by the NMR mechanism called relaxation of the second kind. Calculation time: seconds.

## Physical / mathematical content

- The example constructs damped spin–cavity generators and varies spin–cavity detuning and cavity loss. It extracts the spin-amplitude decay mode from the Liouvillian and converts it to a population relaxation rate, then compares the resonant rate with the exact two-level Purcell expression.

## Numerical / algorithmic content

- In the `zeeman-liouv` formalism without basis approximation, the code diagonalises the Liouvillian across 301 detunings and four loss rates. It also propagates a spin-excitation density operator at selected detunings to compare survival curves.

## Implementation structure

- The spin and cavity isotopes are `{'E','C3'}`; the coupling is `0.35e6`. Detuning spans `2*pi*linspace(-8e6,8e6,301)`, and loss rates are `2*pi*[2e6 4e6 8e6 16e6]`. At the second loss rate, the resonant extracted rate is checked against the exact expression; survival is plotted for detunings `2*pi*[0 2e6 6e6]` over 0–40 μs.
