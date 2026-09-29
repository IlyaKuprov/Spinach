# kernel/evolution.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/evolution.m) | [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=evolution.m)

- Signature: `answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)`

## Purpose and inputs

`evolution` propagates a supplied initial state under the supplied Hamiltonian or Liouvillian. It supports Hilbert-space density matrices and wavefunctions, and Liouville-space state vectors; the active formalism in `spin_system.bas.formalism` selects the implementation. In Liouville calculations, `rho` may be one state vector or a horizontal stack of initial state vectors. In Hilbert-space calculations it is the initial density matrix or wavefunction. In Liouville-space propagation, `coil` is the detection state (one column per channel for `multichannel`) and optional `destination` screens destination states. In Hilbert density-matrix propagation, `coil` is an observable operator and `destination` is ignored. In `zeeman-wavef`, `coil` is a reference wavefunction: `observable` and `multichannel` return overlap trajectories, not operator expectation values. Wavefunction `trajectory` rejects a stack of initial columns, and `total` is unavailable for unitary wavefunction evolution.

The source defines `timestep` as the duration of one step in seconds and `nsteps` as the number of steps. The requested total duration is therefore `timestep*nsteps`. This function has no frequency-offset or field-grid input and does not convert an offset axis to Hz: any offsets or frequencies must already be represented in `L` according to the caller's Hamiltonian/Liouvillian convention. For a manually assembled Liouville generator, the source documents `L=H+1i*R+1i*K`, where `H`, `R`, and `K` are the Hamiltonian-commutation, relaxation, and kinetics superoperators. For Hilbert-space propagation, `L` is the Hamiltonian matrix.

## Propagation direction and rule

Propagation proceeds forward from `rho`. In the ordinary matrix path, the code obtains a one-step propagator using `propagator(spin_system,L,timestep)` and repeatedly left-applies it to the current state, `rho_loc=P*rho_loc`, for the requested steps. Polyadic generators are forwarded to `krylov` with the same time-step and output arguments. The source delegates construction of the propagator to those helpers; it does not define a separate frequency grid or reverse-time propagation here. For Liouville calculations, automatic trajectory-level restriction and subspace splitting can reduce the working space; these are controlled by the Spinach system settings. In a non-Krylov Liouville path the code may choose an internal optimal substep count and set its step to `(timestep*nsteps)/nsteps_opt`; this preserves the requested total duration but does not preserve the caller's individual step spacing. The Krylov path receives the supplied `timestep` and `nsteps` directly.

## Output forms and dimensions

The accepted `output` strings are `final`, `trajectory`, `total`, `refocus`, `observable`, and `multichannel`.

- `final` returns the final state, or a horizontal stack of final state vectors.
- `trajectory` returns the state trajectory. The source rejects this option for a stack of Liouville state vectors; representation depends on the active formalism.
- `observable` returns an observable time trace as a vector for one initial state or a matrix for a stack of initial states.
- `multichannel` returns several observable traces. With a stack of initial states, the documented shape is channels-by-time-by-states; with one initial state it is a channel-by-time matrix.
- `total` integrates the observable trace from the start of the simulation to infinity and requires relaxation. The source rejects it for unitary wavefunction evolution and for polyadic generators.
- `refocus` propagates successive input vectors for zero, one, two, and subsequent steps, matching the described indirect-dimension refocusing use.

The routine returns states or observable values, not a separate time vector; the caller's `timestep` supplies the spacing. The input contract checks that `L`, `coil`, `timestep`, and `nsteps` are numeric, that `rho` is numeric or a cell array, and that `output` is one of the listed strings.

## Parallel context

The source reports Hilbert-space parallel tests through a 128-core configuration (16 nodes, 8 cores each) and notes trajectory parallelisation did not appear beneficial because of inter-thread communication ([10.1063/1.3679656](https://doi.org/10.1063/1.3679656)).
