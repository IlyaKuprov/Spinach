# kernel/trajan.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/trajan.m>

## Purpose

Trajectory analysis function. Plots the time dependence of the density matrix norm, partitioned into user-specified property classes. The trajectory would usually come out of an `evolution.m` run from a given starting point under a given Liouvillian.

## Behaviour

- Syntax: `trajan(spin_system,traj,property,time_axis)`.
- Validates inputs via an internal `grumble` function: the spin system formalism must be `sphten-liouv`; the trajectory must be numeric with a number of rows matching the basis dimension; the property string must be one of the five supported options; if supplied, `time_axis` must be a real row vector with one element per trajectory column.
- Unless the property is `level_populations`, the unit state is projected out of the trajectory before analysis (`traj=traj-(unit*unit')*traj`), so the unit state population is ignored.
- Supported properties:
  - `correlation_order` — time dependence of total populations of one-spin, two-spin, three-spin, etc. correlations. Correlation order of each state is `sum(logical(spin_system.bas.basis),2)`; zero-spin order is eliminated; legend entries are `N-spin`; label `correlation order amplitude`; title `correlation orders`.
  - `coherence_order` — time dependence of different orders of coherence, where coherence order is the sum of projection quantum numbers in the spherical tensor representation of each state (`[~,M]=lin2lm(spin_system.bas.basis)` then `sum(M,2)`). Legend entries are the order in TeX dollar math; label `coherence order amplitude`; title `coherence orders`.
  - `total_each_spin` — time dependence of total state space population involving each individual spin in any way (all local populations and coherences of the spin as well as all of its correlations to third-party spins). Subspace mask is `spin_system.bas.basis(:,n)~=0`; label `density touching each spin`; title `spin populations`.
  - `local_each_spin` — time dependence of the population of the subspace of states local to each individual spin, with no correlations to other spins. Subspace mask is `spin_system.bas.basis(:,n)~=0` combined with `sum(spin_system.bas.basis,2)==spin_system.bas.basis(:,n)`; label `density local to each spin`; title `spin populations`.
  - `level_populations` — populations of the Zeeman energy levels. The trajectory is transformed with `sphten2zeeman(spin_system)*traj` and divided by `prod(spin_system.comp.mults)`; the number of levels is `sqrt(size(traj,1))`; populations are the real parts of the diagonal of each reshaped column; label `energy level populations`; title `level populations`.
- For `total_each_spin` and `local_each_spin`, legend entries use user-specified labels from `spin_system.comp.labels` when available, otherwise isotope names with spin numbers, e.g. `isotope (n)`.
- For all properties except `level_populations`, each result row is the Euclidean norm over the subspace: `sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1))`.
- Y-axis extents: for norm-based properties, `[-0.05*max_val 1.05*max_val]`; for `level_populations`, `[min_val-0.05*(max_val-min_val), max_val+0.05*(max_val-min_val)]`, with a flat trajectory padded by `max(0.05*abs(max_val),1e-6)` so the limits differ.
- Plotting: if a non-empty `time_axis` is supplied, results are plotted against it; otherwise against trajectory point index with x-label `trajectory point`. Line colours are set deterministically as `hsv2rgb([n/numel(p) 0.75 0.75])`. A legend is drawn only if none exists on the current axes (to avoid expensive redraws), using `klegend` with `Location` `Best` and `AutoUpdate` `off`. Axis labels, limits (`xlim tight`), title and grid are applied via `kylabel`, `ylim`, `ktitle` and `kgrid`.
- An unknown property string raises `unknown property.`; the validation helper raises specific errors for wrong formalism, non-numeric trajectory, dimension mismatch, unknown property, or invalid `time_axis`.
- The function is only applicable to trajectories recorded in the sphten-liouv formalism.

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
