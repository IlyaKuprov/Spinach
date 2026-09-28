# examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m

- Signature: `bubble_pulse_acquire()`

## Purpose

Simulates a PNL (partially negative line) spectrum of ortho-deuterium in the presence of a parahydrogenation catalyst. Bubbling is followed by a 45-degree pulse. The source notes that a paper link will follow in due course.

## Physical / mathematical content

- The spin system contains four `2H` nuclei. Experimental chemical shifts are `[4.55, 4.55, -13.5, -16.5]`; J-couplings are 12.0 between spins 1 and 2 and 0.24 between spins 3 and 4.
- NQI tensors for spins 3 and 4 come from a DFT calculation and are rotated from the Gaussian `abc` principal-axis frame into the standard-orientation frame used for the coordinates. The coordinates for spins 1 and 2 are unspecified (D2); those for spins 3 and 4 are `[-1.962, 0.573, -0.576]` and `[-0.175, 1.399, -1.630]`.
- Chemical kinetics use parts `{[1 2], [3 4]}`, rates `[-1 5000; 1 -5000]`, and initial concentrations `[1 0]`. The magnetic field is 7.05.
- Relaxation combines `redfield` and `t1_t2`, with zero equilibrium, secular terms retained, correlation times `{1e-12, 400e-12}`, R1 rates `{0.04 0.04 0 0}`, and R2 rates `{8.00 8.00 0 0}`.

## Numerical / algorithmic content

- The simulation uses the `sphten-liouv` formalism without approximation. It constructs the Hamiltonian, relaxation superoperator, and free-kinetics generator, then adds a bubbling pump based on the deuterium-pair states `S` and `Q{1}` through `Q{5}`. The pump rate `1e-1` is marked in the source as a guess requiring a proper rate.
- Starting from `unit_state`, the system evolves under the Hamiltonian, relaxation, and pumped kinetics for 7 seconds. A pulse-acquire calculation then applies a `pi/4` pulse about `Ly` and detects `L+` on spins 1 and 2.
- Acquisition uses offset `209.6554`, sweep `60`, `256` points, `1024`-point zero filling, ppm axis units, and an inverted axis. The FID receives exponential apodisation with parameter `6`; its FFT is shifted and normalised by the maximum absolute spectral amplitude.

## Implementation structure

- Creates the spin system and basis, runs relaxation analysis, and obtains the deuterium-pair spin states.
- Builds free and bubbling kinetics, evolves the initial state, and passes the resulting state to `liquid` with `hp_acquire`.
- Plots the real spectrum with the y-axis label `NMR intensity, a.u.` and x-axis limits `[4.40 4.70]`.
