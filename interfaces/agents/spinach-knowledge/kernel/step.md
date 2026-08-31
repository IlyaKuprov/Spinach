# kernel/step.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/step.m`
- Signature: `rho=step(spin_system,L,rho,time_step)`
- Total lines: 419

## Purpose

Propagation step function. Computes the action by a matrix exponential without compuing that exponential. Supports one-, two-, and three-point product quadratures. Syntax: rho=step(spin_system,L,rho,time_step)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `comm_series()`, `reordered_taylor()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 52-53: If L{1} is a function handle, route to iserstep; implemented by `if iscell(L)&&(numel(L)==3)&&isa(L{1},'function_handle')`.
- Lines 55-56: Call state-dependent evolution generator solver; implemented by `rho=iserstep(spin_system,L,rho,time_step); return`.
- Lines 60-61: Check consistency; implemented by `grumble(L,rho,time_step)`.
- Lines 63-64: Zero time step shortcut; implemented by `if time_step==0, return; end`.
- Lines 66-67: Do we want to run on GPU?; implemented by `want_gpu=ismember('gpu',spin_system.sys.enable)`.
- Lines 69-70: Is the generator already on GPU?; implemented by `if isnumeric(L)`.
- Lines 82-83: Is the state already on GPU?; implemented by `if isnumeric(rho)`.
- Lines 95-98: Is it expm(A)*x (wavefunctions or state vectors)?; implemented by `expm_times_vec=ismember(spin_system.bas.formalism,{'sphten-liouv', 'zeeman-liouv', 'zeeman-wavef'})`.
- Lines 99-100: expm(A)*rho*expm(-A); implemented by `if ~expm_times_vec`.
- Lines 102-103: Process Lie generators; implemented by `if iscell(L)&&(numel(L)==2)`.
- Lines 105-106: Two-point Lie quadrature; implemented by `L=isergen(L{1},[],L{2},time_step)`.
- Lines 110-111: Three-point Lie quadrature; implemented by `L=isergen(L{1},L{2},L{3},time_step)`.
- Lines 115-116: Fast bypass for small density matrices; implemented by `if size(L,1)<spin_system.tols.small_matrix`.
- Lines 118-119: Use Matlab's expm; implemented by `P=expm(-1i*L*time_step)`.
- Lines 121-122: A cell array; implemented by `if iscell(rho)`.
- Lines 124-125: Parallel over cells; implemented by `parfor n=1:numel(rho)`.
- Lines 131-132: A martrix; implemented by `rho=P*rho*P'`.
- Lines 136-137: Shortcut; implemented by `return`.

### Control flow inferred from the code

- Line 53: conditional branch on `iscell(L)&&(numel(L)==3)&&isa(L{1},'function_handle')`.
- Line 64: conditional branch on `time_step==0, return; end`.
- Line 70: conditional branch on `isnumeric(L)`.
- Line 71: conditional branch on `want_gpu&&(~isa(L,'gpuArray'))`.
- Line 75: `for` loop over `n=1:numel(L)`.
- Line 76: conditional branch on `want_gpu&&(~isa(L{n},'gpuArray'))`.
- Line 83: conditional branch on `isnumeric(rho)`.
- Line 84: conditional branch on `want_gpu&&(~isa(rho,'gpuArray'))`.
- Line 88: `for` loop over `n=1:numel(rho)`.
- Line 89: conditional branch on `want_gpu&&(~isa(rho{n},'gpuArray'))`.
- Line 100: conditional branch on `~expm_times_vec`.
- Line 103: conditional branch on `iscell(L)&&(numel(L)==2)`.
- Line 116: conditional branch on `size(L,1)<spin_system.tols.small_matrix`.
- Line 122: conditional branch on `iscell(rho)`.

### Key state/data transformations

- Lines 56: computes `rho` using `rho=iserstep(spin_system,L,rho,time_step); return`.
- Lines 67: computes `want_gpu` using `want_gpu=ismember('gpu',spin_system.sys.enable)`.
- Lines 72: computes `L` using `L=gpuArray(L)`.
- Lines 77: computes `L{n}` using `L{n}=gpuArray(L{n})`.
- Lines 90: computes `rho{n}` using `rho{n}=gpuArray(rho{n})`.
- Lines 96-98: computes `expm_times_vec` using `expm_times_vec=ismember(spin_system.bas.formalism,{'sphten-liouv', 'zeeman-liouv', 'zeeman-wavef'})`.
- Lines 119: computes `P` using `P=expm(-1i*L*time_step)`.
- Lines 142: computes `norm_gen` using `norm_gen=cheap_norm(L)*abs(time_step); nsteps=ceil(norm_gen/2)`.
- Lines 162: computes `return_numeric` using `return_numeric=true(); rho={rho}`.
- Lines 172-173: computes `parfor_makes_sense` using `parfor_makes_sense=(size(rho{1},1)>256)&&(numel(rho)>32)&& (poolsize>0)&&(~want_gpu)`.
- Lines 182: computes `rho{k}` using `rho{k}=comm_series(spin_system,L,rho{k},time_step,nsteps)`.
- Lines 210: computes `scaling` using `scaling=max(abs(rho),[],'all')`.
- Lines 223: computes `norm_mat` using `norm_mat=cheap_norm(L)*abs(time_step)`.
- Lines 227: computes `nsteps` using `nsteps=ceil(norm_mat/2)`.
- Lines 253: computes `rho(:,m)` using `rho(:,m)=reordered_taylor(L,rho(:,m),time_step,nsteps)`.

### Local helper functions

- Line 272: `comm_series()` — `function rho=comm_series(spin_system,L,rho,time_step,nsteps)`. Convergence tolerance
  - Representative operation: `tol=eps('double')`.
  - Representative operation: `rho=clean_up(spin_system,rho,tol)`.
- Line 328: `reordered_taylor()` — `function rho=reordered_taylor(L,rho,t,nsteps)`. Loop over substeps
  - Representative operation: `for n=1:nsteps`.
  - Representative operation: `next_term=rho; k=1; tol=eps('double')`.
- Line 378: `grumble()` — `function grumble(L,rho,time_step)`.
  - Representative operation: `if (~isnumeric(time_step))||(~isscalar(time_step))`.
  - Representative operation: `error('time_step must be a scalar.')`.

## Parameters / inputs

- L -Liouvillian or Hamiltonian to be used for
- propagation; centre point piecewise-constant
- rule if one matrix is supplied, piecewise-
- linear rule if two matrices {left, right}
- are supplied, piecewise-quadratic if three
- matrices {left, midpoint, right} are given.
- If L is assembled manually from Hamiltonian
- commutation superoperator H, relaxation
- superoperator R, and kinetics superoperator
- K, use L=H+1i*R+1i*K
- State-dependent evolution generators are
- supported: if L{1} is a function handle (see
- iserstep.m documentation), L{2} is current
- time, and L{3} is the method (see iserstep.m
- documentation), the problem is routed to a
- an appropriate Lie group solver.
- rho -state vector or density matrix
- time_step -length of the time step to take

## Outputs

- rho -state vector or density matrix
- Note: we initially had a faithful implementation of the Krylov process
- here -subspace, orthogonalisation, projection, etc., but in all
- our testing it was much inferior to the reordered Taylor process
- that is currently implemented below.
- Note: the peculiar sequence of algebraic operations in the code below
- is designed to minimise the memory footprint in large cases.

## Implementation structure

- Propagation step function. Computes the action by a matrix exponential
- without compuing that exponential. Supports one-, two-, and three-point
- product quadratures. Syntax:
- rho=step(spin_system,L,rho,time_step)
- L - Liouvillian or Hamiltonian to be used for
- propagation; centre point piecewise-constant
- rule if one matrix is supplied, piecewise-
- linear rule if two matrices {left, right}
- are supplied, piecewise-quadratic if three
- matrices {left, midpoint, right} are given.
- If L is assembled manually from Hamiltonian
- commutation superoperator H, relaxation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `iscell()`, `iserstep()`, `grumble()`, `ismember()`, `gpuArray()`, `isergen()`, `cheap_norm()`, `report()`, `num2str()`, `evolution()`, `true()`, `false()`, `comm_series()`, `issparse()`, `cellfun()`, `rho()`.
