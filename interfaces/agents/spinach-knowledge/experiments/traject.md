# experiments/traject.m

- Signature: `traj=traject(spin_system,parameters,H,R,K)`

## Purpose

Records the forward evolution of the supplied initial state as a trajectory over the requested number of points.

## Physical / mathematical content

The function combines the supplied Hamiltonian, relaxation, and kinetics matrices as `L = H + 1i*R + 1i*K`. If requested, analytical decoupling is applied to both this Liouvillian and the initial state before propagation.

## Numerical / algorithmic content

The routine checks that `H`, `R`, and `K` are numeric matrices of the same size; that `parameters.sweep` is a positive real scalar and `parameters.npoints` a positive integer; and that the initial state and decoupling list are supplied. Decoupling labels must be character strings for isotopes present in the system, and analytical decoupling requires the `sphten-liouv` formalism. The evolution uses a time step of `1/parameters.sweep` and requests `parameters.npoints-1` propagation steps in trajectory mode.

## Parameters / inputs

- `parameters.sweep` - sweep width, Hz
- `parameters.npoints` - number of points in the trajectory
- `parameters.rho0` - initial state
- `parameters.decouple` - spins to decouple, e.g. {'15N','13C'}
- `H` - Hamiltonian matrix, received from context function
- `R` - relaxation superoperator, received from context function
- `K` - kinetics superoperator, received from context function

## Outputs

- `traj` - system trajectory, a bookshelf stack of state vectors

## Implementation structure

After consistency checks, the function composes `L = H + 1i*R + 1i*K`, applies the requested decoupling to `L` and `parameters.rho0`, then calls `evolution` with the reciprocal sweep width, the initial state, and trajectory output mode.
