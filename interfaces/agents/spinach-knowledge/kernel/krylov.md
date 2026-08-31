# kernel/krylov.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/krylov.m`
- Signature: `answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)`
- Total lines: 234

## Purpose

Krylov propagation function. Avoids matrix exponentiation, but can be slow. Should be used when the Liouvillian exponential does not fit in- to the system memory, but the Liouvillian itself does. Syntax: answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 70-71: Check consistency; implemented by `grumble(spin_system,L,coil,rho,timestep,nsteps,output)`.
- Lines 73-74: Upload data to GPU or optimise layout; implemented by `if ismember('gpu',spin_system.sys.enable)`.
- Lines 81-82: Start feedback timer; implemented by `feedback=tic()`.
- Lines 84-85: Decide the output type; implemented by `switch output`.
- Lines 89-90: Take the step, ignoring user settings; implemented by `rho=step(spin_system,L,rho,timestep*nsteps)`.
- Lines 92-93: Assign the answer; implemented by `answer=gather(rho)`.
- Lines 97-98: Preallocate the answer and set the starting point; implemented by `answer=zeros([size(rho,1) (nsteps+1)],'like',1i)`.
- Lines 101-102: Loop over steps; implemented by `for n=1:nsteps`.
- Lines 104-105: Take next step; implemented by `rho=step(spin_system,L,rho,timestep)`.
- Lines 107-108: Assign the answer; implemented by `answer(:,n+1)=gather(rho)`.
- Lines 110-111: Inform the user; implemented by `if (n==nsteps)||(toc(feedback)>1)`.
- Lines 121-122: Number of steps is fixed; implemented by `nsteps=size(rho,2)`.
- Lines 124-125: Loop over steps; implemented by `for n=2:nsteps`.
- Lines 127-128: Take next step; implemented by `rho(:,n:end)=step(spin_system,L,rho(:,n:end),timestep)`.
- Lines 144-145: Preallocate the answer; implemented by `answer=zeros([(nsteps+1) size(rho,2)],'like',1i)`.
- Lines 147-148: Set the initial point; implemented by `answer(1,:)=gather(coil'*rho)`.
- Lines 156-157: Assign the answer; implemented by `answer(n+1,:)=gather(coil'*rho)`.
- Lines 170-171: Preallocate the answer; implemented by `answer=zeros([size(coil,2) (nsteps+1)],'like',1i)`.

### Control flow inferred from the code

- Line 74: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 85: dispatches on `output`; cases `'final'`, `'trajectory'`, `'refocus'`, `'observable'`, `'multichannel'`.
- Line 102: `for` loop over `n=1:nsteps`.
- Line 111: conditional branch on `(n==nsteps)||(toc(feedback)>1)`.
- Line 125: `for` loop over `n=2:nsteps`.
- Line 131: conditional branch on `(n==nsteps)||(toc(feedback)>1)`.
- Line 151: `for` loop over `n=1:nsteps`.
- Line 160: conditional branch on `(n==nsteps)||(toc(feedback)>1)`.
- Line 177: `for` loop over `n=1:nsteps`.
- Line 186: conditional branch on `(n==nsteps)||(toc(feedback)>1)`.

### Key state/data transformations

- Lines 75: computes `L` using `L=gpuArray(L); rho=gpuArray(full(rho))`.
- Lines 76: computes `coil` using `coil=gpuArray(coil); location='GPU'`.
- Lines 78: computes `location` using `location='CPU'; rho=full(rho)`.
- Lines 82: computes `feedback` using `feedback=tic()`.
- Lines 90: computes `rho` using `rho=step(spin_system,L,rho,timestep*nsteps)`.
- Lines 93: computes `answer` using `answer=gather(rho)`.
- Lines 99: computes `answer(:,1)` using `answer(:,1)=gather(rho)`.
- Lines 108: computes `answer(:,n+1)` using `answer(:,n+1)=gather(rho)`.
- Lines 122: computes `nsteps` using `nsteps=size(rho,2)`.
- Lines 128: computes `rho(:,n:end)` using `rho(:,n:end)=step(spin_system,L,rho(:,n:end),timestep)`.
- Lines 148: computes `answer(1,:)` using `answer(1,:)=gather(coil'*rho)`.
- Lines 157: computes `answer(n+1,:)` using `answer(n+1,:)=gather(coil'*rho)`.

### Local helper functions

- Line 204: `grumble()` — `function grumble(spin_system,L,coil,rho,timestep,nsteps,output)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv', 'zeeman-wavef'})`.
  - Representative operation: `'zeeman-wavef'})`.

## Parameters / inputs

- L -the Liouvillian to be used during evolution
- rho -the initial state vector or a horizontal stack thereof
- output -a string giving the type of evolution that is required
- 'final' -returns the final state vector or a horizontal
- stack thereof.
- 'trajectory' -returns the stack of state vectors giving
- the trajectory of the system starting from
- rho with the user-specified number of steps
- and step length.
- 'total' -returns the integral of the observable trace
- from the simulation start to infinity. This
- option requires the presence of relaxation.
- 'refocus' -evolves the first vector for zero steps,
- second vector for one step, third vector for
- two steps, etc., consistent with the second
- stage of evolution in the indirect dimension
- after a refocusing pulse.
- 'observable' -returns the time dynamics of an observable
- as a vector (if starting from a single ini-
- tial state) or a matrix (if starting from a
- stack of initial states).
- 'multichannel' -returns the time dynamics of several
- observables as rows of a matrix. Note
- that destination state screening may be
- less efficient when there are multiple
- destinations to screen against.
- coil -the detection state, used when 'observable' is specified as
- the output option. If 'multichannel' is selected, the coil
- should contain multiple columns corresponding to individual
- observable vectors.

## Outputs

- answer -a vector or a matrix, depending on the options set during
- the call.
- Note: this function does not support the zeeman-hilb formalism; in
- zeeman-wavef, L is the Hamiltonian matrix, rho is a wavefunc-
- tion or a horizontal stack thereof, coil is a reference wave-
- function, and observables are overlap trajectories.
- Note: we initially had a faithful implementation of the Krylov process
- here -subspace, orthogonalisation, projection, etc., but in all
- our testing it was much inferior to the reordered Taylor process
- that is currently implemented below.

## Implementation structure

- Krylov propagation function. Avoids matrix exponentiation, but can be
- slow. Should be used when the Liouvillian exponential does not fit in-
- to the system memory, but the Liouvillian itself does. Syntax:
- answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)
- L -the Liouvillian to be used during evolution
- rho -the initial state vector or a horizontal stack thereof
- output -a string giving the type of evolution that is required
- 'final' -returns the final state vector or a horizontal
- stack thereof.
- 'trajectory' -returns the stack of state vectors giving
- the trajectory of the system starting from
- rho with the user-specified number of steps

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `gpuArray()`, `tic()`, `step()`, `gather()`, `answer()`, `toc()`, `report()`, `num2str()`, `rho()`, `iscell()`, `ischar()`.
