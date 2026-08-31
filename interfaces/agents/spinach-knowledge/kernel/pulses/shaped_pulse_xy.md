# kernel/pulses/shaped_pulse_xy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/shaped_pulse_xy.m`
- Signature: `[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...`
- Total lines: 409

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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 67-68: Check consistency; implemented by `grumble(drift,controls,amplitudes,slice_durs,rho,method)`.
- Lines 70-71: Start feedback timer; implemented by `feedback=tic()`.
- Lines 73-74: Get trajectory started; implemented by `if nargout>1`.
- Lines 76-77: Preallocate as a cell array; implemented by `traj=cell(1,numel(slice_durs)+1)`.
- Lines 79-81: First element is initial condition; implemented by `traj{1}=rho`.
- Lines 85-86: Get the propagator started; implemented by `if nargout>2, P=speye(size(drift)); end`.
- Lines 88-89: Move appropriate things to GPU; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 98-99: Decide propagation method; implemented by `switch method(1:4)`.
- Lines 101-102: Krylov method; implemented by `case 'expv'`.
- Lines 104-105: Loop over pulse slices; implemented by `for n=1:numel(slice_durs)`.
- Lines 107-108: Decide the quadrature; implemented by `switch method(6:8)`.
- Lines 110-111: Piecewise-constant; implemented by `case 'pwc'`.
- Lines 113-114: Grab the drift; implemented by `slice_oper=drift`.
- Lines 116-117: Add controls; implemented by `for k=1:numel(controls)`.
- Lines 121-122: Apply the evolution slice; implemented by `rho=step(spin_system,slice_oper,rho,slice_durs(n))`.
- Lines 124-125: Update the propagator; implemented by `if nargout>2`.
- Lines 130-131: 2nd order Lie; implemented by `case 'pwl'`.
- Lines 133-134: Grab left and right edge drifts; implemented by `slice_oper_l=drift; slice_oper_r=drift`.

### Control flow inferred from the code

- Line 74: conditional branch on `nargout>1`.
- Line 86: conditional branch on `nargout>2, P=speye(size(drift)); end`.
- Line 89: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 91: `for` loop over `n=1:numel(controls)`.
- Line 95: conditional branch on `nargout>2, P=gpuArray(P); end`.
- Line 99: dispatches on `method(1:4)`; cases `'expv'`, `'pwc'`, `'pwl'`.
- Line 105: `for` loop over `n=1:numel(slice_durs)`.
- Line 108: dispatches on `method(6:8)`; cases `'pwc'`, `'pwl'`.
- Line 117: `for` loop over `k=1:numel(controls)`.
- Line 125: conditional branch on `nargout>2`.
- Line 137: `for` loop over `k=1:numel(controls)`.
- Line 146: conditional branch on `nargout>2`.
- Line 159: conditional branch on `nargout>1, traj{n+1}=rho; end`.
- Line 162: conditional branch on `(n==numel(slice_durs))||(toc(feedback)>1)`.

### Key state/data transformations

- Lines 71: computes `feedback` using `feedback=tic()`.
- Lines 77: computes `traj` using `traj=cell(1,numel(slice_durs)+1)`.
- Lines 81: computes `traj{1}` using `traj{1}=rho`.
- Lines 90: computes `drift` using `drift=gpuArray(drift)`.
- Lines 92: computes `controls{n}` using `controls{n}=gpuArray(controls{n})`.
- Lines 94: computes `rho` using `rho=gpuArray(rho)`.
- Lines 114: computes `slice_oper` using `slice_oper=drift`.
- Lines 126: computes `P` using `P=step(spin_system,slice_oper,P,slice_durs(n))`.
- Lines 134: computes `slice_oper_l` using `slice_oper_l=drift; slice_oper_r=drift`.
- Lines 139: computes `slice_oper_r` using `slice_oper_r=slice_oper_r+amplitudes{k}(n+1)*controls{k}`.
- Lines 217: computes `PS` using `PS=propagator(spin_system,slice_oper,slice_durs(n))`.
- Lines 340: computes `traj{n}` using `traj{n}=gather(traj{n})`.

### Local helper functions

- Line 349: `grumble()` — `function grumble(drift,controls,amplitudes,slice_durs,rho,method)`.
  - Representative operation: `if (~isnumeric(drift))||(size(drift,1)~=size(drift,2))`.
  - Representative operation: `error('drift operator must be a square matrix.')`.

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
