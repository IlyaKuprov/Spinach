# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_conversion.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_conversion.m`
- Signature: `mq_conversion()`
- Total lines: 110

## Purpose

Optimal control design of the multiple-quantum conversion pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the conversion pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is three rotor periods (240 us) long in 0.5 us slices, the controls are Cartesian, and the 100 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is the Hermitian combination of the +MQ and -MQ coherences between the m=+3/2 and m=-3/2 levels (3Q) or between the m=+5/2 and m=-5/2 levels (5Q), and the target is the population difference across the central transition, as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation; the waveforms supplied in this folder reached fidelities of 0.91 (3Q) and 0.86 (5Q) after 500 iterations.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control design of the multiple-quantum conversion pulse
- of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
- nach, the conversion pulse optimisation from
- A single 27Al nucleus with the quadrupolar coupling and the shi-
- elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
- ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
- magnet. The quadrupolar interaction is taken to second order in
- the rotating frame, and the powder average runs over 200 crystal-
- lite orientations at 80 initial rotor phases each. The pulse is
- three rotor periods (240 us) long in 0.5 us slices, the controls
- are Cartesian, and the 100 kHz amplitude ceiling is enforced by
- a spillout penalty followed by clipping. The initial state is the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `assume()`, `mqmas_drifts()`, `operator()`, `randi()`, `optimcon()`, `fmaxnewton()`, `cartesian2polar()`, `polar2cartesian()`, `grape_xy()`.
