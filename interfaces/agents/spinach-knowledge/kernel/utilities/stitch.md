# kernel/utilities/stitch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/stitch.m`
- Signature: `fid=stitch(spin_system,L,rho_stack,coil_stack,...`
- Total lines: 200

## Purpose

Stitching function for bidirectionally propagated 3D NMR pulse sequences. Propagate your initial condition forward to some mid- point, propagate your detection state backward to the same mid- point and use this function to obtain the 3D free induction de- cay (http://dx.doi.org/10.1016/j.jmr.2014.04.002) at the price of two 2D simulations. Syntax: fid=stitch(spin_system,L,rho_stack,coil_stack,... mtp_oper,mtp_time,t1

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Set default time directions; implemented by `if ~exist('tdir','var'), tdir='+-'; end`.
- Lines 52-53: Check consistency; implemented by `grumble(L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3)`.
- Lines 55-56: Preallocate the fid; implemented by `fid=zeros(t3.nsteps,t2.nsteps,t1.nsteps,'like',1i)`.
- Lines 58-59: Compute time step propagators; implemented by `P=propagator(spin_system,L,t2.timestep/2); Pct=P'`.
- Lines 61-62: Compute midpoint event propagator; implemented by `Pm=speye(size(L))`.
- Lines 68-69: Decide the device and style; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 71-72: Upload the problem to GPU; implemented by `P=gpuArray(P)`.
- Lines 78-79: Inform the user; implemented by `report(spin_system,'stitching will be done on GPU.')`.
- Lines 83-84: Inform the user; implemented by `report(spin_system,'stitching will be done on CPU.')`.
- Lines 88-89: Stitch the trajectories; implemented by `for k=1:t2.nsteps`.
- Lines 91-93: Inform the user; implemented by `report(spin_system,['stitching forward and backward trajectories, step ' num2str(k) '/' num2str(t2.nsteps) ' '])`.
- Lines 95-96: Scalar product with the midpoint propagator; implemented by `fid(:,k,:)=gather(coil_stack'*(Pm*rho_stack))`.
- Lines 98-99: Time directions; implemented by `switch tdir`.
- Lines 103-104: Forward time evolution for rho; implemented by `rho_stack=P*rho_stack`.
- Lines 106-107: Forward time evolution for coil; implemented by `coil_stack=P*coil_stack`.
- Lines 114-115: Backward time evolution for coil; implemented by `coil_stack=Pct*coil_stack`.
- Lines 119-120: Backward time evolution for rho; implemented by `rho_stack=Pct*rho_stack`.
- Lines 135-136: Complain and bomb out; implemented by `error('invalid time direction specification in tdir')`.

### Control flow inferred from the code

- Line 50: conditional branch on `~exist('tdir','var'), tdir='+-'; end`.
- Line 63: `for` loop over `n=1:numel(mec_oper)`.
- Line 69: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 89: `for` loop over `k=1:t2.nsteps`.
- Line 99: dispatches on `tdir`; cases `'++'`, `'+-'`, `'-+'`, `'--'`.

### Key state/data transformations

- Lines 56: computes `fid` using `fid=zeros(t3.nsteps,t2.nsteps,t1.nsteps,'like',1i)`.
- Lines 59: computes `P` using `P=propagator(spin_system,L,t2.timestep/2); Pct=P'`.
- Lines 62: computes `Pm` using `Pm=speye(size(L))`.
- Lines 74: computes `Pct` using `Pct=gpuArray(Pct)`.
- Lines 75: computes `rho_stack` using `rho_stack=gpuArray(rho_stack)`.
- Lines 76: computes `coil_stack` using `coil_stack=gpuArray(coil_stack)`.
- Lines 96: computes `fid(:,k,:)` using `fid(:,k,:)=gather(coil_stack'*(Pm*rho_stack))`.

### Local helper functions

- Line 145: `grumble()` — `function grumble(L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3)`.
  - Representative operation: `if (~isnumeric(L))||(~ismatrix(L))||(size(L,1)~=size(L,2))`.
  - Representative operation: `error('L must be a square matrix.')`.

## Parameters / inputs

- L -spin system Liouvillian
- rho_stack -state vector stack from the forward part of
- the simulation
- coil_stack -coil vector stack from the backward part of
- the sumulation
- mec_oper -cell array of operators in the midpoint event
- chain, e.g. {Lx,L,Sy}
- mec_time -cell array of durations of each event at the
- midpoint of the t2 evolution period
- t1.nsteps -number of time steps in t1
- t2.nsteps -number of time steps in t2
- t2.timestep -duration of each time step in t2
- t3.nsteps -number of time steps in t3
- tdir -time direction for state and coil propagation,
- the default is '+-'

## Outputs

- fid -three-dimensional free induction decay

## Implementation structure

- Stitching function for bidirectionally propagated 3D NMR pulse
- sequences. Propagate your initial condition forward to some mid-
- point, propagate your detection state backward to the same mid-
- point and use this function to obtain the 3D free induction de-
- cay (http://dx.doi.org/10.1016/j.jmr.2014.04.002) at the price
- of two 2D simulations. Syntax:
- fid=stitch(spin_system,L,rho_stack,coil_stack,...
- mtp_oper,mtp_time,t1,t2,t3)
- L -spin system Liouvillian
- rho_stack -state vector stack from the forward part of
- the simulation
- coil_stack -coil vector stack from the backward part of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `propagator()`, `speye()`, `clean_up()`, `ismember()`, `gpuArray()`, `report()`, `num2str()`, `fid()`, `gather()`, `ismatrix()`, `isfield()`, `isscalar()`, `iscell()`, `all()`.
