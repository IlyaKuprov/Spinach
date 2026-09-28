# kernel/trajsimil.m

- Signature: `score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)`

## Purpose

Compares two state-space trajectories at corresponding time points and returns a similarity-score vector. For further information, see [the JMR article](http://dx.doi.org/10.1016/j.jmr.2013.02.012).

## Parameters / inputs

- `spin_system` - Spinach system whose basis and formalism define the trajectory states.
- `trajectory_1`, `trajectory_2` - equal-size numeric trajectory matrices, with basis states in rows and time points in columns; their row count must match the basis size.
- `scorefcn` - one of `'RSP'`, `'RDN'`, or a grouped form: `'SG-RSP'`, `'SG-RDN'`, `'BSG-RSP'`, or `'BSG-RDN'`.
  - `RSP` computes the column-wise scalar product with `dot`.
  - `RDN` computes `1 - norm(trajectory_1(:,k)-trajectory_2(:,k),2)/2` for each time point `k`.
  - The `SG-` prefix groups `T(l,m)` and `T(l,-m)` states as equivalent; `BSG-` groups all non-identity states of each spin. For grouped scores, corresponding grouped contributions are formed by summing coefficient absolute-squares and taking the square root.

## Output

- `score` - one similarity score per time point.

## Behavior

The function checks that the trajectories have equal size and match the basis row count. Ungrouped scoring requires `sphten-liouv` or `zeeman-liouv` formalism; the `SG-` and `BSG-` options require `sphten-liouv`.
