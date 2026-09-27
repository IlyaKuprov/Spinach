# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m

- Signature: `ct_selective()`

## Purpose

Optimal control design of the central transition selective pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the soft pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is 50 us long in 0.5 us slices, the controls are Cartesian, and the 10 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is the population difference across the central transition, and the target is the single-quantum coherence of the central transition, as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation; the waveform supplied in this folder reached a fidelity of 0.69 after 500 iterations, against the maximum of 1/sqrt(2) for this initial and target state pair.

## Implementation

- Builds the 27Al Zeeman-Hilbert spin system and obtains MAS drift Hamiltonians over the 200-orientation two-angle grid and 80 initial rotor phases (160 rotor ticks).
- The normalized initial state is the central-transition population difference; the target is central-transition single-quantum coherence.
- Optimizes Cartesian controls in 100 slices of 0.5 µs (50 µs total), with a 10 kHz amplitude ceiling, SNSA spillout penalty, L-BFGS, and a 500-iteration limit. The amplitude is clipped and the fidelity recomputed afterward.
- Saves the waveform in rad/s and slice durations to ct_pulse.mat.
