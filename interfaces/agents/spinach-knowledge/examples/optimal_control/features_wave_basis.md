# examples/optimal_control/features_wave_basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_wave_basis.m`
- Signature: `features_wave_basis()`
- Total lines: 116

## Purpose

An illustration of basis set coefficient optimisation for a pul- se optimised as a linear combination of user-specified vectors. The optimisation is done using LBFGS GRAPE algorithm. In a set of 100 equspaced signals, the central 60 spins are set up for maximum excitation; there are no constraints on the dyna- mics of the 20 spins on either side of the interval. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- An illustration of basis set coefficient optimisation for a pul-
- se optimised as a linear combination of user-specified vectors.
- The optimisation is done using LBFGS GRAPE algorithm.
- In a set of 100 equspaced signals, the central 60 spins are set
- up for maximum excitation; there are no constraints on the dyna-
- mics of the 20 spins on either side of the interval.
- Calculation time: minutes.
- Set the magnetic field
- 100 non-interacting spins at equal intervals
- within the plus/minus 160 ppm range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `wave_basis()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`.
