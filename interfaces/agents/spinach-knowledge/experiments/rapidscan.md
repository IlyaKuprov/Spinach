# experiments/rapidscan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/rapidscan.m`
- Signature: `[b_axis,spectrum]=rapidscan(spin_system,parameters)`
- Total lines: 120

## Purpose

Time-domain rapid field scan ESR experiment, Eatons style. Syntax: [b_axis,spectrum]=rapidscan(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `carrier()`, `rotframe()`, `state()`, `spectrum()`, `step()`, `waveform()`, `isfield()`, `isscalar()`.
