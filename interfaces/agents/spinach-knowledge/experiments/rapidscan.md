# experiments/rapidscan.m

- Signature: `[b_axis,spectrum]=rapidscan(spin_system,parameters)`

## Purpose

Simulates a time-domain rapid-scan ESR experiment in the electron rotating frame. It combines the Zeeman and coupling terms with the microwave drive and relaxation, then propagates the equilibrium state as the magnetic field is swept.

## Numerical / algorithmic content

The routine starts from isotropic thermal equilibrium, constructs the `L+` detection state, and advances the state for `parameters.nsteps` time steps using `step`. At each step it records the `L+` observable amplitude and uses the corresponding field offset in the propagator. The returned field axis is the sweep waveform shifted by the centre field `spin_system.inter.magnet`.

## Parameters / inputs

- `parameters.mw_pwr` — microwave power in rad/s.
- `parameters.sweep` — two ascending magnetic-field sweep offsets in Tesla, relative to the centre field in `spin_system.inter.magnet`.
- `parameters.nsteps` — number of magnetic-field steps.
- `parameters.timestep` — duration of each time step in seconds.

## Outputs

- `b_axis` — magnetic-field axis in Tesla.
- `spectrum` — `L+` observable amplitude at each magnetic field.

Call this experiment directly, without a context.
