# kernel/evolution.m

- Signature: `answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)`

## Purpose

Time evolution function. Performs all types of time propagation with automatic trajectory level state space restriction. Syntax: answer=evolution(spin_system,L,coil,rho,timestep,... nsteps,output,destination)

## Physical / mathematical content

- Supports Hilbert-space density-matrix, Liouville-space state-vector, and wavefunction propagation, with trajectory-level state-space restriction.

## Numerical / algorithmic content

- Polyadic generators are forwarded to `krylov()`. In the Liouville-space `final` path, smaller subspaces use an exponential propagator and large subspaces use Krylov propagation; the `krylov` enable/disable settings can affect this choice.
- Hilbert-space final-state and observable calculations support parallel execution; the source reports tests through 128 cores and says parallel trajectory calculation did not appear beneficial because of inter-thread communication.

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

- `answer` is a vector, matrix, channels-by-time-by-states array, or cell array of matrices, depending on the selected output and formalism.
- The source reports Hilbert-space parallel tests through 128 cores and notes trajectory parallelization did not appear beneficial because of inter-thread communication ([10.1063/1.3679656](https://doi.org/10.1063/1.3679656)).
