# examples/optimal_control/state_transfer_m2s.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_m2s.m`
- Signature: `state_transfer_m2s()`
- Total lines: 147

## Purpose

A transfer of coherence from longitudinal magnetization into a two-spin singlet state in allyl pyruvate with a distribution of B1 powers and transmitter offsets. XX and YY components of the singlet dephase rapidly in this system, and are therefore drop- ped from the target state specification. Calculation time: many hours.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A transfer of coherence from longitudinal magnetization into a
- two-spin singlet state in allyl pyruvate with a distribution of
- B1 powers and transmitter offsets. XX and YY components of the
- singlet dephase rapidly in this system, and are therefore drop-
- ped from the target state specification.
- Calculation time: many hours.
- Get the spin system
- Kill the methyl group
- Magnetic field (500.13 MHz)
- Formalism and basis set
- Run Spinach housekeeping
- Set up and normalise the initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `allyl_pyruvate()`, `create()`, `basis()`, `state()`, `idxof()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `cumsum()`, `fmaxnewton()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`.
