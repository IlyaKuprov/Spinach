# kernel/trajan.m

- Signature: `trajan(spin_system,traj,property,time_axis)`

## Purpose

Plots trajectory analysis results. For further information, see [the JMR article](http://dx.doi.org/10.1016/j.jmr.2013.02.012).

## Parameters / inputs

- `spin_system` - Spinach system and basis description.
- `traj` - numeric trajectory matrix with one basis state per row and one time point per column.
- `property` - one of:
  - `'correlation_order'` - norm of the trajectory in each nonzero spin-correlation-order subspace.
  - `'coherence_order'` - norm in each coherence-order subspace; the order is the sum of projection quantum numbers in the spherical-tensor representation.
  - `'total_each_spin'` - norm of all states involving each spin, including its local states and correlations with other spins.
  - `'local_each_spin'` - norm of the states local to each spin, excluding correlations with other spins.
  - `'level_populations'` - real Zeeman-basis diagonal populations, normalized by the product of spin multiplicities.
- `time_axis` - optional numeric real row vector with one time value per trajectory column. If omitted or empty, trajectory-point indices are used.

## Output

Writes the selected plot into the current figure. Requires trajectories recorded in `sphten-liouv` formalism.

## Behavior

For all properties except `'level_populations'`, the unit-state component is projected out before analysis. The function then plots the norm or populations specified by `property`; each property is checked against the supported options above.
