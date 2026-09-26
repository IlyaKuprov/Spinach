# experiments/rapidscan.m

- Signature: `[b_axis,spectrum]=rapidscan(spin_system,parameters)`

## Purpose

Time-domain rapid field scan ESR experiment, Eatons style. Syntax: [b_axis,spectrum]=rapidscan(spin_system,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.mw_pwr -microwave power, rad/s
- parameters.sweep -magnetic field sweep extents,
- arouind the centre field specified
- in sys.magnet, a two-element vector
- in Tesla
- parameters.nsteps -number of steps in the magnetic
- field sweep
- parameters.timestep -duration of each time step, seconds

## Outputs

- b_axis -magnetic field axis, Tesla
- spectrum -L+ observable amplitude at each
- magnetic field
- Note: this experiment should be called directly without a context.

## Implementation structure

- Time-domain rapid field scan ESR experiment, Eatons style. Syntax:
- [b_axis,spectrum]=rapidscan(spin_system,parameters)
- parameters.mw_pwr -microwave power, rad/s
- parameters.sweep -magnetic field sweep extents,
- arouind the centre field specified
- in sys.magnet, a two-element vector
- in Tesla
- parameters.nsteps -number of steps in the magnetic
- field sweep
- parameters.timestep -duration of each time step, seconds
- b_axis -magnetic field axis, Tesla
- spectrum -L+ observable amplitude at each
