# experiments/eqmag.m

- Signature: `magn=eqmag(spin_system,parameters)`

## Purpose

Computes the molar magnetization vector at the thermal equilibrium at the temperature specified in inter.temperature and magnetic field spe- cified in sys.magnet (assumed to be along the Z-axis), averaged over system orientations using the spherical grid specified. Syntax: magn=eqmag(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.grid -spherical grid for averaging

## Outputs

- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.

## Implementation structure

- Computes the molar magnetization vector at the thermal equilibrium at
- the temperature specified in inter.temperature and magnetic field spe-
- cified in sys.magnet (assumed to be along the Z-axis), averaged over
- system orientations using the spherical grid specified. Syntax:
- magn=eqmag(spin_system,parameters)
- parameters.grid -spherical grid for averaging
- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.
- Check consistency
