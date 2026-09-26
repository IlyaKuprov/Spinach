# experiments/fieldscan_magn.m

- Signature: `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`

## Purpose

Z magnetization of the sample as a function of magnetic field in a finite-speed magnetic field sweep experiment. Syntax: [fields,z_magn]=fieldscan_magn(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
