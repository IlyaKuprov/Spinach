# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m`
- Signature: `ct_selective()`
- Total lines: 105

## Purpose

Optimal control design of the central transition selective pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the soft pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is 50 us long in 0.5 us slices, the controls are Cartesian, and the 10 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is the population difference across the central transition, and the target is the single-quantum coherence of the central transition, as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation; the waveform supplied in this folder reached a fidelity of 0.69 after 500 iterations, against the maximum of 1/sqrt(2) for this initial and target state pair.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control design of the central transition selective pulse
- of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
- nach, the soft pulse optimisation from
- A single 27Al nucleus with the quadrupolar coupling and the shi-
- elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
- ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
- magnet. The quadrupolar interaction is taken to second order in
- the rotating frame, and the powder average runs over 200 crystal-
- lite orientations at 80 initial rotor phases each. The pulse is
- 50 us long in 0.5 us slices, the controls are Cartesian, and the
- 10 kHz amplitude ceiling is enforced by a spillout penalty follo-
- wed by clipping. The initial state is the population difference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `assume()`, `mqmas_drifts()`, `operator()`, `randi()`, `optimcon()`, `fmaxnewton()`, `cartesian2polar()`, `polar2cartesian()`, `grape_xy()`.
