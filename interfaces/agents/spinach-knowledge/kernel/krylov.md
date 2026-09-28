# kernel/krylov.m

- Signature: `answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)`

## Purpose

Propagates one or more state vectors without forming the full propagator. The routine uses a reordered Taylor process rather than a Krylov-subspace/Arnoldi iteration.

## Physical / mathematical content

For `sphten-liouv` and other supported Liouville-space use, `L` is the evolution generator and `rho` contains the initial state vector or vectors. In `zeeman-wavef`, `L` is the Hamiltonian, `rho` contains wavefunctions, and observables are overlap trajectories. The `zeeman-hilb` formalism is not supported.

## Numerical / algorithmic content

Each time step is evaluated by the Taylor-based `step` routine; GPU arrays are used when GPU execution is enabled in `spin_system.sys.enable`. The method avoids storing a full matrix exponential and may be slow.

## Parameters / inputs

- `spin_system` - Spinach system descriptor and execution settings.
- `L` - Liouvillian or, in `zeeman-wavef`, Hamiltonian used for evolution.
- `coil` - detection state for `observable`; for `multichannel`, columns are the individual observable vectors. In `zeeman-wavef`, it is a reference wavefunction.
- `rho` - initial state vector or horizontal stack of initial states; in `zeeman-wavef`, wavefunction(s).
- `timestep` - time step.
- `nsteps` - number of time steps.
- `output` - output mode: `final` returns the state after `nsteps`; `trajectory` returns the initial and subsequent states; `refocus` evolves successive input vectors for zero, one, two, and further steps; `observable` returns the coil signal versus time; `multichannel` returns signals for multiple coil vectors.

## Outputs

- `answer` - final state, state trajectory, or observable trajectory, with dimensions depending on `output` and whether `rho` or `coil` contains multiple columns. `multichannel` returns one channel per coil vector; destination-state screening can be less efficient with multiple destinations.

## Implementation structure

The routine validates its inputs, moves data to the GPU when enabled, dispatches on `output`, and advances the requested state or observable trajectories through repeated calls to `step`.
