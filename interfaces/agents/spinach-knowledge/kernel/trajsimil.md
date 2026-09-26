# kernel/trajsimil.m

- Signature: `score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)`

## Purpose

Computes trajectory similarity scores. Returns a function representing "similarity" of the two state space trajectories at different points in time. See http://dx.doi.org/10.1016/j.jmr.2013.02.012 for further infor- mation. Syntax: score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Parameters / inputs

- trajectory_1,2 -spin system trajectories, supplied as nstates
- x nsteps matrices.
- scorefcn -similarity scoring method; possibilities are:
- 'RSP' -running scalar product. Computes
- scalar products between the cor-
- responding vectors of the trajec-
- tories.
- 'RDN' -running difference norm. The two
- trajectories are subtracted and
- difference 2-norms returned.
- 'SG-' -prefix that turns on state grou-
- ping. T(l,m) and T(l,-m) states
- of each spin (standalone or in
- direct products with other ope-
- rators)will be considered equva-
- lent.
- 'BSG-' -prefix that turns on broad state
- grouping. All states of a given
- spin (standalone or in direct
- products with other operators)
- will be considered equivalent.
- The possible combinations are: 'RSP','RDN',
- 'SG-RSP','SG-RDN','BSG-RSP','BSG-RDN'.
- State grouping consists in summing the absolute squares of the coeffi-
- cients to be grouped and taking the square root. The trajectories would
- usually come out of the evolution.m or krylov.m run from a given star-
- ting point under a given Liouvillian.
- Output:
- score -the similarity score vector, one element
- per time slice
- Note: SG and BSG options require sphten-liouv formalism.

## Implementation structure

- Computes trajectory similarity scores. Returns a function representing
- "similarity" of the two state space trajectories at different points in
- time. See http://dx.doi.org/10.1016/j.jmr.2013.02.012 for further infor-
- mation. Syntax:
- score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)
- trajectory_1,2 -spin system trajectories, supplied as nstates
- x nsteps matrices.
- scorefcn -similarity scoring method; possibilities are:
- 'RSP' -running scalar product. Computes
- scalar products between the cor-
- responding vectors of the trajec-
- tories.
