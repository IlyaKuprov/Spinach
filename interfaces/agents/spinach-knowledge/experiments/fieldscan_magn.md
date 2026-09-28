# experiments/fieldscan_magn.m

- Signature: `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`

## Purpose

Z magnetization of the sample as a function of magnetic field in a finite-speed magnetic field sweep experiment. Syntax: [fields,z_magn]=fieldscan_magn(spin_system,parameters)

## Physical / mathematical content


- At the specified orientation, constructs the z magnetic-moment operator from the rotated g tensors and spin operators.
- Initializes the density operator at the first field and evaluates magnetization with `hdot(rho,mz)`.

## Numerical / algorithmic content


- Creates a linearly spaced magnetic-field grid and adjusts the Zeeman Hamiltonian across the sweep.
- If `parameters.nstates` is supplied, the requested eigenstates are used to form an active-space projection; the magnetization observable is the real part of `hdot(rho,mz)`.

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


- Validates the `zeeman-hilb` setup and sweep inputs, then constructs the field grid and orientation-dependent magnetic-moment operator.
- Initializes the density operator, optionally forms an active-space projection, and evaluates the magnetization observable in the acquisition loop.
