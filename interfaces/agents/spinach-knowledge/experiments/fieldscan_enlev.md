# experiments/fieldscan_enlev.m

- Signature: `fieldscan_enlev(spin_system,parameters)`

## Purpose

Plots a user-specified number of the lowest energy levels of the system as a function of the applied magnetic field. The energies are obtained using the Arnoldi method. Syntax: fieldscan_enlev(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Parameters / inputs

- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.nstates -number of lowest energy states
- to solve for

## Implementation structure

- Plots a user-specified number of the lowest energy levels of
- the system as a function of the applied magnetic field. The
- energies are obtained using the Arnoldi method. Syntax:
- fieldscan_enlev(spin_system,parameters)
- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.nstates -number of lowest energy states
