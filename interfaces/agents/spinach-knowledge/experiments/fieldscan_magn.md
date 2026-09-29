# experiments/fieldscan_magn.m

- MATLAB implementation: [experiments/fieldscan_magn.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/fieldscan_magn.m)

Source: https://spindynamics.org/wiki/index.php?title=fieldscan_magn.m

Signature: [fields,z_magn]=fieldscan_magn(spin_system,parameters).

This simulates the sample's z magnetisation during a finite-speed magnetic-field sweep. For the fixed requested orientation, it constructs lab-frame Zeeman and coupling Hamiltonians and a z magnetic-moment observable from the rotated g tensors and spin operators. The initial density operator is thermal equilibrium at the starting field, normalised to unit trace. If parameters.nstates is supplied, the calculation is projected into the corresponding low-energy subspace.

The field grid is linear between the requested endpoints. With npoints samples and total sweep_time, the time increment is dt=sweep_time/(npoints-1). At each sample the code records real(hdot(rho,mz)) and propagates the state with the propagator for the instantaneous field-dependent Hamiltonian. The source does not provide a physical unit for the returned magnetisation values, so retain them as the routine's simulated signal rather than assigning a unit.

Required inputs:

- parameters.fields: two ascending endpoints in tesla.
- parameters.npoints: number of field samples; at least two are needed for the source's npoints-1 time-step denominator.
- parameters.orientation: three Euler angles [alp bet gam] in radians.
- parameters.sweep_time: total sweep duration in seconds; the source requires it to be positive.
- spin_system: must use zeeman-hilb formalism. The optional parameters.nstates is a positive integer specifying the low-energy active-space size.

Return values: fields is the sampled magnetic-field axis in tesla and z_magn is the corresponding simulated magnetisation signal. This is not, by itself, a DNP/hyperpolarisation protocol or an imaging acquisition: the source contains no RF irradiation, polarisation-transfer step, spatial encoding, or explicit relaxation/kinetics input during the sweep.

Limit: the model is the finite-rate coherent field sweep from a thermal starting state at one fixed orientation, with optional low-energy projection. It is not a claim about measured magnetisation or experimental enhancement.
