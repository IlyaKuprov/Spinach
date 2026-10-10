# kernel/trajsimil.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/trajsimil.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/trajsimil.m)

## Purpose

Computes trajectory similarity scores. Returns one numeric similarity score per time column for two state-space trajectories.

## Behaviour

- Syntax: `score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)`.
- Consistency is enforced first: the function is only available for Liouville space formalisms (`sphten-liouv` or `zeeman-liouv`); `scorefcn` must be a character string; the two trajectory matrices must have matching dimensions; each trajectory's row count must equal the basis set dimension; both trajectories must be numeric arrays of doubles; and `scorefcn` must be one of `RSP`, `RDN`, `SG-RSP`, `SG-RDN`, `BSG-RSP`, `BSG-RDN`.
- If `scorefcn` starts with `SG-` or `BSG`, state grouping is run before scoring:
  - `SG-`: all T(l,-m) states are renamed into T(l,m) states via `lin2lm` and `lm2lin` with `abs(M)`.
  - `BSG`: all non-identity states are renamed into Lz (entries of the state list not equal to 0 are set to 2).
  - Grouping is performed separately inside each substance block, so identical local descriptors (including units) of different substances are never combined. Unique grouped states are found with `unique(...,'rows')`, and for each group the trajectory rows are combined by root-sum-square: `sqrt(sum(abs(...).^2,1))` over the coefficients belonging to that group.
  - A progress message reports collapsing equivalent subspaces and, when finished, how many states were collected into how many groups.
  - State grouping (SG and BSG) is only available for the `sphten-liouv` formalism.
- After any grouping, the score function is computed per time slice:
  - `RSP` (running scalar product): `score(n)=dot(trajectory_1(:,n),trajectory_2(:,n))` for each column.
  - `RDN` (running difference norm): `score(n)=1-norm(trajectory_1(:,n)-trajectory_2(:,n),2)/2` for each column.
  - Any other value raises the error `unknown similarity score function.`
- The trajectories would usually come out of `evolution.m` or `krylov.m` run from a given starting point under a given Liouvillian.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object; its compiled dimension (`spin_system.bas.offsets(end)`) must match the trajectory row dimension, and its formalism must be a Liouville space formalism.
- `trajectory_1`, `trajectory_2` — spin system trajectories, supplied as `nstates x nsteps` matrices; both must be arrays of doubles with identical dimensions.
- `scorefcn` — similarity scoring method:
  - `'RSP'` — running scalar product; computes scalar products between the corresponding vectors of the trajectories.
  - `'RDN'` — running difference norm; the two trajectories are subtracted and difference 2-norms returned.
  - `'SG-'` — prefix that turns on state grouping; T(l,m) and T(l,-m) states of each spin (standalone or in direct products with other operators) are considered equivalent.
  - `'BSG-'` — prefix that turns on broad state grouping; all states of a given spin (standalone or in direct products with other operators) are considered equivalent.
  - Possible combinations: `'RSP'`, `'RDN'`, `'SG-RSP'`, `'SG-RDN'`, `'BSG-RSP'`, `'BSG-RDN'`.

**Outputs**

- `score` — the similarity score vector, one element per time slice (a `1 x nsteps` row vector).

## References

- Article: [http://dx.doi.org/10.1016/j.jmr.2013.02.012](http://dx.doi.org/10.1016/j.jmr.2013.02.012)
- Wiki: [https://spindynamics.org/wiki/index.php?title=trajsimil.m](https://spindynamics.org/wiki/index.php?title=trajsimil.m)
