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
- observables as rows of a matrix (if
- starting from a single initial state)
- or as a channels-by-time-by-states
- array (if starting from a stack of
- initial states). Note that destination
- state screening may be less efficient
- when there are multiple destinations
- to screen against.
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

- answer -a vector, a matrix, a channels-by-time-by-states array,
- or a cell array of matrices, depending on the options
- set during the call
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
