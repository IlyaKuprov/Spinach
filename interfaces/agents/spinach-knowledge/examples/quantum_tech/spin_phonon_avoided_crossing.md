# examples/quantum_tech/spin_phonon_avoided_crossing.m

- Signature: `spin_phonon_avoided_crossing()`

## Purpose

Shows the avoided crossing between an electron-spin transition and a quantised phonon mode in a resonant spin–phonon exchange model. The phonon is requested with the `V#` particle syntax. Calculation time: seconds.

## Physical / mathematical content

- The one-excitation spin–phonon doublet is diagonalised while the spin–phonon detuning is swept. Resonant exchange splits the dressed levels; the minimum gap is checked against the coupling-derived value.

## Numerical / algorithmic content

- The `zeeman-hilb` model uses no basis approximation. The code isolates the one-quantum manifold, diagonalises the projected Hamiltonian at 121 detunings, and plots the two dressed energies.

## Implementation structure

- The system is `{'E','V3'}`, with the phonon at zero rotating-frame frequency and exchange coupling `4e6`. The detuning grid is `2*pi*linspace(-20e6,20e6,121)`; the resonant gap is validated against `2*g/(2*pi*1e6)` in MHz.
