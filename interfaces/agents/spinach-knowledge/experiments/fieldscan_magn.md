# experiments/fieldscan_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/fieldscan_magn.m`
- Signature: `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`
- Total lines: 162

## Purpose

Z magnetization of the sample as a function of magnetic field in a finite-speed magnetic field sweep experiment. Syntax: [fields,z_magn]=fieldscan_magn(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
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
- parameters.sweep_time -sweep time, seconds
- parameters.nstates -(optional) number of lowest energy
- states to use for the effective
- Hamiltonian in the time domain

## Outputs

- fields -magnetic fields in Tesla at each point in time
- z_magn -total sample magnetisation at each point in time
- Note: this function requires Hilbert space formalism.

## Implementation structure

- Z magnetization of the sample as a function of magnetic field in a
- finite-speed magnetic field sweep experiment. Syntax:
- [fields,z_magn]=fieldscan_magn(spin_system,parameters)
- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.sweep_time -sweep time, seconds
- parameters.nstates -(optional) number of lowest energy

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `equilibrium()`, `orientation()`, `isfield()`, `z_magn()`, `hdot()`, `propagator()`, `fields()`, `strcmp()`, `isrow()`.
