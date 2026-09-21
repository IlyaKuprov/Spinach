# examples/optimal_control/distortions/distortions_figure_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_3.m`
- Signature: `distortions_figure_3()`
- Total lines: 225

## Purpose

Figure 3 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Figure 3 from the paper by Rasulov and Kuprov:
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping
- Set up spin states
- Get the control operators
- Get the drift Hamiltonian
- Define control parameters
- Last five slices are dead time

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `guess()`, `fmaxnewton()`, `xy_profile()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `figure()`.
