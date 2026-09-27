# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m

- Signature: `mq_excitation()`

## Purpose

Optimal control design of the multiple-quantum excitation pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the excitation pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is three rotor periods (240 us) long in 0.5 us slices, the controls are Cartesian, and the 100 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is Iz and the target is the Hermitian combination of the +MQ and -MQ coherences between the m=+3/2 and m=-3/2 levels (3Q) or between the m=+5/2 and m=-5/2 levels (5Q), as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation; the waveforms supplied in this folder reached fidelities of 0.72 (3Q) and 0.65 (5Q) after 500 iterations.

## Implementation

- Uses the 27Al MAS drift model on a 200-orientation powder grid, with 80 initial rotor phases and 160 rotor ticks.
- The normalized initial state is Lz; the target is the Hermitian combination of the +MQ and −MQ coherences for the selected 3Q or 5Q order. The script selects 5Q by default.
- Optimizes 480 Cartesian slices of 0.5 µs (240 µs total) with a 100 kHz ceiling, SNSA spillout penalty, L-BFGS, and up to 500 iterations; it clips the amplitude and reevaluates fidelity.
- Saves the waveform and slice durations as mq_exc_<order>q.mat.
