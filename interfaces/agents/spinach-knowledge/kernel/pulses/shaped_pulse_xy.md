# kernel/pulses/shaped_pulse_xy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/shaped_pulse_xy.m`
- Signature: `[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...`
- Total lines: 422

## Purpose

Shaped pulse function using Cartesian coordinates. Applies a user- specified pulse shape on user-specified operators while the rest of the drift Liouvillian continues to affect the spin system. Syntax: [rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,... amplitudes,slice_durs,rho,method)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- drift -the drift Liouvillian, the part of the Liouvillian that
- should continue running in the background. This should
- include the transmitter offset term, if any.
- controls -a cell array of control operators corresponding to each
- channel, this may include operators for spatial degrees
- of freedom, such as gradients and diffusion.
- amplitudes -a cell array of control amplitude vectors in rad/s, one
- vector per control channel; the elements of each vector
- correspond to different time points.
- slice_durs -a vector containing the duration of each pulse slice,
- seconds. For piecewise-constant methods, the number of
- durations should be equal to the nuber of amplitudes.
- For piecewise-linear methods, there should be one ele-
- ment more in the amplitude array.
- rho -initial state vector or a bookshelf matrix thereof
- method -propagation method and product quadrature:
- Krylov algorithm (usually faster for calls with one
- and two outputs):
- 'expv-pwc' -piecewise-constant
- 'expv-pwl' -2nd order Lie quadrature
- Explicit matrix exponentiation (usually faster for
- calls with three outputs):
- 'expm-pwc' -piecewise-constant
- 'expm-pwl' -2nd order Lie quadrature
- Spinach evolution function call (do not choose un-
- less you have a specific good reason):
- 'evol-pwc' -piecewise-constant
- 'evol-pwl' -2nd order Lie quadrature

## Outputs

- rho -state vector for the final state, or a stack thereof
- traj -system trajectory as a [1 x (nsteps+1)] cell array,
- the first point is the initial condition
- P -effective pulse propagator (expensive, best avoided)

In `zeeman-hilb`, all six methods return the one-sided ordered propagator `P`, reusable as `P*rho*P'`; the returned state and trajectory still use two-sided density-matrix evolution. Requesting this third output with `expv-*` or `evol-*` explicitly constructs slice propagators without changing the state-propagation method.

## Implementation structure

- Shaped pulse function using Cartesian coordinates. Applies a user-
- specified pulse shape on user-specified operators while the rest of
- the drift Liouvillian continues to affect the spin system. Syntax:
- [rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...
- amplitudes,slice_durs,rho,method)
- drift -the drift Liouvillian, the part of the Liouvillian that
- should continue running in the background. This should
- include the transmitter offset term, if any.
- controls -a cell array of control operators corresponding to each
- channel, this may include operators for spatial degrees
- of freedom, such as gradients and diffusion.
- amplitudes -a cell array of control amplitude vectors in rad/s, one

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `tic()`, `speye()`, `ismember()`, `gpuArray()`, `method()`, `step()`, `slice_durs()`, `clean_up()`, `toc()`, `report()`, `num2str()`, `isergen()`, `propagator()`, `evolution()`, `gather()`.
