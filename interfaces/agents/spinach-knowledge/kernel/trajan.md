# kernel/trajan.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/trajan.m>

## Purpose

Trajectory analysis function. Plots the time dependence of the density matrix norm, partitioned into user-specified property classes. The trajectory would usually come out of an `evolution.m` run from a given starting point under a given Liouvillian.

## Interpretation of plotted trajectories

For `sphten-liouv` state-vector columns, each substance unit coordinate is removed independently before analysing correlation or coherence content; `level_populations` is the exception. The non-population plots show norms of selected coefficient subspaces, not a time-resolved probability distribution.

- `correlation_order` groups basis states by the number of spins carrying a nontrivial factor; `coherence_order` groups them by the sum of spherical-tensor projection quantum numbers. Each trace is the norm of the corresponding coefficient subspace.
- `total_each_spin` includes every state involving that spin, including correlations with other spins. `local_each_spin` retains only states local to that spin. Because total-spin groups can overlap, their plotted norms should not be summed as disjoint populations.
- `level_populations` transforms to the Zeeman representation, divides each substance block by its own Hilbert dimension and plots its real diagonal, concatenating levels in substance order; the unit-state contribution is retained.

Descriptor selections use local spin columns and the compiled offsets; spin-resolved plots retain global spin numbering.

The optional `time_axis` locates the trajectory columns on the plot. The function writes to the current figure and returns no numerical array.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object; its basis must be in `sphten-liouv` formalism.
- `traj` — a stack of state vectors of any length; the number of rows must match the number of states in the basis.
- `property` — one of `correlation_order`, `coherence_order`, `total_each_spin`, `local_each_spin`, `level_populations`.
- `time_axis` — optional user-specified time axis; a row vector of time positions of each state vector in the trajectory array.

**Outputs**

- This function writes into the current figure; no variables are returned.

## References

- Article behind the method: <http://dx.doi.org/10.1016/j.jmr.2013.02.012>
- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=trajan.m>
