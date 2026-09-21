# experiments/fieldscan_enlev.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/fieldscan_enlev.m`
- Signature: `fieldscan_enlev(spin_system,parameters)`
- Total lines: 113

## Purpose

Plots a user-specified number of the lowest energy levels of the system as a function of the applied magnetic field. The energies are obtained using the Arnoldi method. Syntax: fieldscan_enlev(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `orientation()`, `hz2icm()`, `kfigure()`, `kxlabel()`, `kylabel()`, `strcmp()`, `isfield()`, `isrow()`.
