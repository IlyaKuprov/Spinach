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

## Outputs

- answer -a vector, a matrix, or a channels-by-time-by-states
- array, depending on the options set during the call.
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
