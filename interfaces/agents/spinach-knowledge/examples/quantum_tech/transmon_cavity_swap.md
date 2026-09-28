# examples/quantum_tech/transmon_cavity_swap.m

- Signature: `transmon_cavity_swap()`

## Purpose

Vacuum Rabi swap between a transmon and a microwave cavity mode, both represented by truncated bosonic Weyl algebras. This is the circuit-QED Jaynes-Cummings limit of Blais et al., Rev. Mod. Phys. 93, 025005 (2021). Calculation time: seconds.

## Model and parameters

At zero magnet field, the rotating-frame transmon (T3) and cavity (C3) frequencies are both zero; the transmon anharmonicity is -250 MHz and their exchange coupling is 20 MHz. The calculation uses the Zeeman-Hilbert formalism with no approximation and starts with the transmon in BL2 and the cavity in BL1.

## Calculation

The cavity-context device trajectory uses `sweep=2e9` and `npoints=301`. The code projects transmon and cavity excitation populations, checks that exchange is visible (cavity population reaches at least 0.95 and transmon population falls to at most 0.05), and verifies active-doublet population conservation to 1e-6. It plots the populations over 0–150 ns. The model is cited to Blais et al., Rev. Mod. Phys. 93, 025005 (2021).
