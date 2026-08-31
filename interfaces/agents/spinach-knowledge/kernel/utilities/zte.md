# kernel/utilities/zte.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/zte.m`
- Signature: `projector=zte(spin_system,L,rho,nstates)`
- Total lines: 174

## Purpose

Zero track elimination function. Inspects the first few steps in the system trajectory and drops the states that did not get populated to a user-specified tolerance. Syntax: projector=zte(spin_system,L,rho,nstates)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Validate the input; implemented by `grumble(spin_system,L,rho)`.
- Lines 48-49: Run Zero Track Elimination; implemented by `if ismember('zte',spin_system.sys.disable)`.
- Lines 51-52: Skip if instructed to do so by the user; implemented by `report(spin_system,'WARNING - zero track elimination disabled, basis left unchanged.')`.
- Lines 54-55: Return a unit matrix; implemented by `projector=1`.
- Lines 59-60: Skip if the benefit is likely to be minor; implemented by `report(spin_system,'WARNING - too few zeros in the state vector, basis left unchanged.')`.
- Lines 67-68: Skip if the state vector norm is too small for Krylov procedure; implemented by `report(spin_system,'WARNING - state vector norm below drop tolerance, basis left unchanged.')`.
- Lines 75-76: Get the time step; implemented by `timestep=1/cheap_norm(L)`.
- Lines 78-79: Do not allow infinite time step; implemented by `if isinf(timestep)`.
- Lines 83-84: Report to the user; implemented by `if exist('nstates','var')`.
- Lines 93-94: Preallocate the trajectory; implemented by `trajectory=zeros(numel(rho),spin_system.tols.zte_nsteps,'like',1i)`.
- Lines 96-97: Set the starting point; implemented by `trajectory(:,1)=rho`.
- Lines 100-101: Compute trajectory steps with Krylov technique; implemented by `for n=2:spin_system.tols.zte_nsteps`.
- Lines 103-104: Take a step forward; implemented by `trajectory(:,n)=step(spin_system,L,trajectory(:,n-1),timestep)`.
- Lines 106-107: Analyze the trajectory; implemented by `prev_space_dim=nnz(max(abs(trajectory(:,1:(n-1))),[],2)>spin_system.tols.zte_tol)`.
- Lines 110-112: Inform the user; implemented by `report(spin_system,['evolution step ' num2str(n-1) ', active space dimension ' num2str(curr_space_dim)])`.
- Lines 114-115: Terminate if done early; implemented by `if curr_space_dim==prev_space_dim, break; end`.
- Lines 119-120: Determine which tracks to drop; implemented by `if exist('nstates','var')`.
- Lines 122-123: Determine state amplitudes; implemented by `amplitudes=max(abs(trajectory),[],2)`.

### Control flow inferred from the code

- Line 49: conditional branch on `ismember('zte',spin_system.sys.disable)`.
- Line 79: conditional branch on `isinf(timestep)`.
- Line 84: conditional branch on `exist('nstates','var')`.
- Line 101: `for` loop over `n=2:spin_system.tols.zte_nsteps`.
- Line 115: conditional branch on `curr_space_dim==prev_space_dim, break; end`.
- Line 120: conditional branch on `exist('nstates','var')`.

### Key state/data transformations

- Lines 55: computes `projector` using `projector=1`.
- Lines 76: computes `timestep` using `timestep=1/cheap_norm(L)`.
- Lines 94: computes `trajectory` using `trajectory=zeros(numel(rho),spin_system.tols.zte_nsteps,'like',1i)`.
- Lines 97: computes `trajectory(:,1)` using `trajectory(:,1)=rho`.
- Lines 104: computes `trajectory(:,n)` using `trajectory(:,n)=step(spin_system,L,trajectory(:,n-1),timestep)`.
- Lines 107: computes `prev_space_dim` using `prev_space_dim=nnz(max(abs(trajectory(:,1:(n-1))),[],2)>spin_system.tols.zte_tol)`.
- Lines 108: computes `curr_space_dim` using `curr_space_dim=nnz(max(abs(trajectory),[],2)>spin_system.tols.zte_tol)`.
- Lines 123: computes `amplitudes` using `amplitudes=max(abs(trajectory),[],2)`.
- Lines 126: computes `[~,index]` using `[~,index]=sort(amplitudes,'descend')`.
- Lines 129: computes `zero_track_mask` using `zero_track_mask=true(size(rho))`.
- Lines 130: computes `zero_track_mask(index(1:nstates))` using `zero_track_mask(index(1:nstates))=false()`.

### Local helper functions

- Line 151: `grumble()` — `function grumble(spin_system,L,rho)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})`.
  - Representative operation: `error('zero track elimination is only available for zeeman-liouv and sphten-liouv formalisms.')`.

## Parameters / inputs

- L -the Liouvillian to be used for time
- propagation
- rho -the initial state to be used for
- time propagation
- nstates -if this parameter is specified, only
- nstates most populated states are kept,
- irrespective of the tolerance parameter
- Output:
- projector -projector matrix into the reduced space,
- to be used as follows:
- L_reduced=P'*L*P
- rho_reduced=P'*rho;
- Note: default tolerance may be altered by setting sys.tols.zte_tol
- variable before calling create.m
- Note: further information on how this function works is available
- in IK's JMR paper on the subject
- Note: if tiny interactions or nearly equivalent spins are present,
- it is best to disable zero track elimination by adding 'zte'
- to the sys.disable cell array.

## Implementation structure

- Zero track elimination function. Inspects the first few steps in the
- system trajectory and drops the states that did not get populated to
- a user-specified tolerance. Syntax:
- projector=zte(spin_system,L,rho,nstates)
- L -the Liouvillian to be used for time
- propagation
- rho -the initial state to be used for
- time propagation
- nstates -if this parameter is specified, only
- nstates most populated states are kept,
- irrespective of the tolerance parameter
- Output:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `report()`, `nnz()`, `cheap_norm()`, `isinf()`, `exist()`, `num2str()`, `trajectory()`, `step()`, `true()`, `zero_track_mask()`, `index()`, `false()`, `speye()`, `projector()`.
