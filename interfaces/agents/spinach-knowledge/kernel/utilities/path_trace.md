# kernel/utilities/path_trace.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/path_trace.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/path_trace.m)

## Purpose

Liouvillian path tracing. Treats the user-supplied Liouvillian as the adjacency matrix of a graph, computes the weakly connected subgraphs of that graph and returns a cell array of projectors into independently evolving populated subspaces.

## Behaviour

- Syntax: `projectors=path_trace(spin_system,L,rho)`.
- Input consistency is enforced by an internal `grumble` subfunction: both `L` and `rho` must be numeric, `L` must be square, and if `rho` is non-empty its row count must match the column count of `L`.
- If `'pt'` appears in `spin_system.sys.disable`, path tracing is disabled and a unit projector `{1}` is returned after a warning.
- If `size(L,2)` is smaller than `spin_system.tols.merge_dim`, path tracing is skipped and a unit projector `{1}` is returned.
- The connectivity matrix is built as `G=(abs(L)>spin_system.tols.liouv_zero)`, then symmetrised with its transpose and combined with the identity so isolated states are not lost: `G=or(G,transpose(G)); G=or(G,speye(size(G)))`.
- Weakly connected subgraphs are obtained with `scomponents(G)`; the number of subspaces is `max(member_states)`.
- If `rho` is non-empty, subspace population screening runs in a `parfor` loop using `spin_system.tols.subs_drop`:
  - For the `'sphten-liouv'` and `'zeeman-liouv'` formalisms (Liouville space, state vectors), a subspace is important when `norm(rho.*(member_states==n),1)>tolerance`.
  - For the `'zeeman-hilb'` formalism (Hilbert space, state matrices), a subspace is important when either `norm(rho(member_states==n,:),1)>tolerance` or `norm(rho(:,member_states==n),1)>tolerance`.
  - Any other formalism specification raises the error `'unexpected formalism specification.'`.
- Projectors into significant subspaces are built as sparse matrices of size `size(L,1)` by the subspace dimension, with columns selecting the member state indices.
- Dimension statistics are reported per unique subspace dimension, followed by the total number of kept subspaces and their total dimension.
- Unless `'merge'` is in `spin_system.sys.disable`, small subspaces are merged into batches using `binpack(subspace_dims,spin_system.tols.merge_dim)`; the projectors in each bin are horizontally concatenated into a single projector, and the resulting working subspace dimensions are reported. If merging is disabled, a warning is printed and the unmerged projectors are kept.
- The returned projectors are intended to be used as `L_reduced=P'*L*P` for matrices and `rho_reduced=P'*rho` for state vectors.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object supplying `spin_system.sys.disable`, `spin_system.bas.formalism`, and the tolerances `spin_system.tols.liouv_zero`, `spin_system.tols.subs_drop`, and `spin_system.tols.merge_dim`.
- `L` — Hamiltonian or Liouvillian matrix; must be numeric and square.
- `rho` — the initial state (source state screening) or the detection state (destination state screening); pass `[]` to disable screening. If non-empty, must be numeric with row count consistent with `L`.

Outputs:

- `projectors` — a cell array of projectors into independently evolving populated subspaces, to be used as `L_reduced=P'*L*P` (matrices) and `rho_reduced=P'*rho` (state vectors).

## References

- [http://dx.doi.org/10.1063/1.3398146](http://dx.doi.org/10.1063/1.3398146)
- [http://dx.doi.org/10.1016/j.jmr.2011.03.010](http://dx.doi.org/10.1016/j.jmr.2011.03.010)
- Spin Dynamics Wiki: [https://spindynamics.org/wiki/index.php?title=path_trace.m](https://spindynamics.org/wiki/index.php?title=path_trace.m)
