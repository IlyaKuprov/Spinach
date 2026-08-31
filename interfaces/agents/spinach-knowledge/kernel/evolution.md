# kernel/evolution.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/evolution.m`
- Signature: `answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)`
- Total lines: 1005

## Purpose

Time evolution function. Performs all types of time propagation with automatic trajectory level state space restriction. Syntax: answer=evolution(spin_system,L,coil,rho,timestep,... nsteps,output,destination)

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

- Lines 118-119: Check consistency; implemented by `grumble(L,coil,rho,timestep,nsteps,output)`.
- Lines 121-122: Call Krylov propagation for polyadics; implemented by `if isa(L,'polyadic')`.
- Lines 127-128: Gather state vectors from GPUs; implemented by `if isa(rho,'gpuArray'), rho=gather(rho); end`.
- Lines 130-131: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 135-136: Apply trajectory-level reduction algorithms; implemented by `report(spin_system,'trying to reduce the problem dimension ')`.
- Lines 138-139: Decide the screening algorithm; implemented by `if ismember('dss',spin_system.sys.disable)`.
- Lines 141-142: If DSS is disabled, run forward reduction; implemented by `report(spin_system,'WARNING - destination state screening is disabled.')`.
- Lines 147-148: If destination state is supplied, use it for DSS; implemented by `report(spin_system,'destination state screening using supplied state.')`.
- Lines 153-154: If coil state is supplied, use it for DSS; implemented by `report(spin_system,'destination state screening using coil state.')`.
- Lines 159-160: Default to the usual forward screening; implemented by `report(spin_system,'destination state screening is not applicable.')`.
- Lines 165-166: Count subspaces; implemented by `nsubs=numel(projectors)`.
- Lines 168-169: Parallel strategy; implemented by `if nsubs>1`.
- Lines 175-176: Run the evolution; implemented by `switch output`.
- Lines 180-181: Create arrays of projections; implemented by `L_sub=cell(nsubs,1); rho_sub=cell(nsubs,1)`.
- Lines 189-190: Loop in parallel over independent subspaces; implemented by `parfor (sub=1:nsubs,nworkers)`.
- Lines 192-193: Inform the user; implemented by `report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) ' '])`.
- Lines 195-196: Grab local copies; implemented by `L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub}`.
- Lines 200-203: Decide if Krylov should be used; implemented by `if ((size(L_loc,2)<spin_system.tols.krylov_tol)|| ismember('krylov',spin_system.sys.disable))&& (~ismember('krylov',spin_system.sys.enable))`.

### Control flow inferred from the code

- Line 122: conditional branch on `isa(L,'polyadic')`.
- Line 128: conditional branch on `isa(rho,'gpuArray'), rho=gather(rho); end`.
- Line 131: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv','zeeman-wavef'}`, `'final'`, `'total'`, `'trajectory'`.
- Line 139: conditional branch on `ismember('dss',spin_system.sys.disable)`.
- Line 169: conditional branch on `nsubs>1`.
- Line 176: dispatches on `output`; cases `'final'`, `'total'`, `'trajectory'`.
- Line 184: `for` loop over `sub=1:nsubs`.
- Line 190: `parfor` loop over `(sub=1:nsubs,nworkers)`.
- Line 201: conditional branch on `((size(L_loc,2)<spin_system.tols.krylov_tol)||`.
- Line 220: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 229: `for` loop over `n=1:nsteps_opt`.
- Line 242: `for` loop over `n=1:nsteps_opt`.
- Line 281: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-wavef')`.
- Line 289: `for` loop over `sub=1:nsubs`.

### Key state/data transformations

- Lines 124: computes `answer` using `answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output); return`.
- Lines 143: computes `projectors` using `projectors=reduce(spin_system,L,rho)`.
- Lines 166: computes `nsubs` using `nsubs=numel(projectors)`.
- Lines 170: computes `nworkers` using `nworkers=poolsize`.
- Lines 181: computes `L_sub` using `L_sub=cell(nsubs,1); rho_sub=cell(nsubs,1)`.
- Lines 182: computes `rows` using `rows=cell(nsubs,1); cols=cell(nsubs,1); vals=cell(nsubs,1)`.
- Lines 185: computes `L_sub{sub}` using `L_sub{sub}=projectors{sub}'*L*projectors{sub}`.
- Lines 186: computes `rho_sub{sub}` using `rho_sub{sub}=projectors{sub}'*rho`.
- Lines 196: computes `L_loc` using `L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub}`.
- Lines 198: computes `rho_loc` using `rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol)`.
- Lines 206: computes `nsteps_opt` using `nsteps_opt=ceil(log2(size(L_loc,1)/(size(rho_loc,2)*log(2))))`.
- Lines 208: computes `timestep_loc` using `timestep_loc=total_time/nsteps_opt`.
- Lines 211-213: computes `report(spin_system,['dim(L)` using `report(spin_system,['dim(L)=' num2str(size(L_loc,1)) ', rho stack size ' num2str(size(rho_loc,2)) ', ' num2str(nsteps_opt) ' substeps in optimal stepping.'])`.
- Lines 217: computes `P` using `P=propagator(spin_system,L_loc,timestep_loc)`.
- Lines 258: computes `[rows{sub},cols{sub},vals{sub}]` using `[rows{sub},cols{sub},vals{sub}]=find(proj_loc*sparse(rho_loc))`.
- Lines 287: computes `rho_sub` using `rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1)`.
- Lines 292: computes `coil_sub{sub}` using `coil_sub{sub}=projectors{sub}'*coil`.
- Lines 308: computes `coil_loc` using `coil_loc=clean_up(spin_system,coil_loc,spin_system.tols.zte_tol)`.

### Local helper functions

- Line 974: `grumble()` — `function grumble(L,coil,rho,timestep,nsteps,output)`.
  - Representative operation: `if ~isnumeric(L)`.
  - Representative operation: `error('Liouvillian must be numeric.')`.

## Parameters / inputs

- For Liouville space calculations:
- L -the Liouvillian to be used during evolution. If L
- is assembled manually from Hamiltonian commutation
- superoperator H, relaxation superoperator R, and
- kinetics superoperator K, use L=H+1i*R+1i*K.
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
- destination -(optional) the state to be used for destination state
- screening.
- For Hilbert space calculations:
- L -Hamiltonian matrix
- coil -observable operator (if any)
- rho -initial density matrix
- timestep -duration of a single time step (seconds)
- nsteps -number of steps to take
- output -a string giving the type of evolution that is required
- 'final' -returns the final density matrix.
- 'trajectory' -returns a cell array of density matrices
- giving the trajectory of the system star-
- ting from rho with the user-specified num-
- ber of steps and step length.
- 'refocus' -evolves the first matrix for zero steps,
- second matrix for one step, third matrix for
- two steps, etc., consistent with the second
- stage of evolution in the indirect dimension
- after a refocusing pulse.
- 'observable' -returns the time dynamics of an observable
- as a vector.
- destination -this argument is ignored.
- For wavefunction calculations the Liouville space call signature
- applies with L the Hamiltonian matrix, rho a wavefunction or a
- horizontal stack thereof, and coil a reference wavefunction: the
- 'observable' and 'multichannel' outputs return overlap trajectories
- of the coil with the evolving wavefunction; expectation values of
- operators require a density matrix formalism. Stacks are supported
- by 'final', 'refocus', 'observable', and 'multichannel'; the
- 'trajectory' output takes a single column, and 'total' is not
- defined for unitary wavefunction evolution.

## Outputs

- answer -a vector, a matrix, or a cell array or matrices,
- depending on the options set during the call
- Calculation of final states and observables in Hilbert space is parallel-
- ized and tested all the way to 128-core (16 nodes, 8 cores each) configu-
- rations. Parallelization of the trajectory calculation does not appear to
- yield any benefits due to large amount of inter-thread communication. See

## Implementation structure

- Time evolution function. Performs all types of time propagation with
- automatic trajectory level state space restriction. Syntax:
- answer=evolution(spin_system,L,coil,rho,timestep,...
- nsteps,output,destination)
- For Liouville space calculations:
- L -the Liouvillian to be used during evolution. If L
- is assembled manually from Hamiltonian commutation
- superoperator H, relaxation superoperator R, and
- kinetics superoperator K, use L=H+1i*R+1i*K.
- rho -the initial state vector or a horizontal stack thereof
- output -a string giving the type of evolution that is required
- 'final' -returns the final state vector or a horizontal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `krylov()`, `gather()`, `ismember()`, `reduce()`, `exist()`, `parfor()`, `num2str()`, `clean_up()`, `log2()`, `dim()`, `propagator()`, `gpuArray()`, `clear()`, `cell2mat()`.
