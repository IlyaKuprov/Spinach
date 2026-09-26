# kernel/trajan.m

- Signature: `trajan(spin_system,traj,property,time_axis)`

## Purpose

Trajectory analysis function. Plots the time dependence of the densi- ty matrix norm, partitioned into user-specified property classes. See for further information. Syntax: trajan(spin_system,traj,property,time_axis)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- traj -a stack of state vectors of any length. The
- number of rows in the trajectory array must
- match the number of states in the basis.
- property -if set to 'correlation_order', returns the
- time dependence of the total populations of
- one-spin, two-spin, three-spin, etc. corre-
- lations in the system.
- if set to 'coherence_order', returns the ti-
- me dependence of different orders of coheren-
- ce in the system, where a coherence order is
- defined as the sum of projection quantum num-
- bers in the spherical tensor representation
- of each state.
- if set to 'total_each_spin', returns the ti-
- me dependence of total state space populati-
- on that involves each individual spin in the
- system in any way (all local populations and
- coherences of the spin as well as all of its
- correlations to all third party spins).
- if set to 'local_each_spin', returns the ti-
- me dependence of the population of the sub-
- space of states that are local to each indi-
- vidual spin and do not involve any correla-
- tions to other spins in the system.
- if set to 'level_populations', returns the
- populations of the Zeeman energy levels.
- time_axis -(optional) user specified time axis, a row
- vector of time positions of each state vec-
- tor inthe trajectory array.
- The trajectory would usually come out of the evolution.m run from a
- given starting point under a given Liouvillian.
- Output:
- this function writes into the current figure
- Note: this function is only applicable to the trajectories recorded
- in sphten-liouv formalism.
- Note: unit state population is ignored.

## Implementation structure

- Trajectory analysis function. Plots the time dependence of the densi-
- ty matrix norm, partitioned into user-specified property classes. See
- for further information. Syntax:
- trajan(spin_system,traj,property,time_axis)
- traj -a stack of state vectors of any length. The
- number of rows in the trajectory array must
- match the number of states in the basis.
- property -if set to 'correlation_order', returns the
- time dependence of the total populations of
- one-spin, two-spin, three-spin, etc. corre-
- lations in the system.
- if set to 'coherence_order', returns the ti-
