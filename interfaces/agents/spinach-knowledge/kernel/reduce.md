# kernel/reduce.m

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/reduce.m
- Wiki: https://spindynamics.org/wiki/index.php?title=reduce.m
- Signature: `projectors=reduce(spin_system,L,rho)`

## Purpose

Returns projectors into independently evolving reduced subspaces, selected using the supplied Liouvillian `L` and state `rho`. This is state-space reduction; this function does not construct a relaxation generator or define relaxation-term units.

The local projector cells in `bas.sym_fact(n)` are embedded at `bas.offsets(n)` before screening. Each projector therefore acts on one substance; no tensor products between substances are formed. These spin-only projectors are used only when the supplied generator has the compiled spin dimension and no reaction records are declared. Reaction maps need not preserve the substance-local irreps; reaction-bearing systems therefore use full-generator zero-track elimination and path tracing, rather than deleting inter-irrep chemical couplings. Enlarged spatial-spin generators instead proceed without spin symmetry factorisation; zero-track elimination and path tracing retain their usual enabled behaviour.

Permutation symmetry separates irreducible sectors; zero-track elimination removes coordinates unoccupied by the trajectory; path tracing separates disconnected subspaces. The available reductions depend on the formalism and the system's enable and disable flags.

## Numerical / algorithmic content

Horizontal stacks of wavefunctions or Liouville states retain their actual complex columns during symmetry screening. The matrix 1-norm of a projected stack is its largest column 1-norm, so occupancy is tested without averaging columns or discarding their phases. Liouville stacks continue through ZTE and path tracing without disabling useful reduction. ZTE propagates state columns in bounded multi-column batches, combining row-wise amplitude maxima; column-wise normalisation prevents weak states from disappearing under shared Taylor convergence.

Hilbert-space density matrices in cell arrays are screened one by one in their original phase; any matrix that occupies an irrep keeps that sector, without averaging or magnitude proxies. Input validation precedes reduction.

## Syntax

```matlab
projectors=reduce(spin_system,L,rho)
```

## Parameters / inputs

- `spin_system`: Spinach system, including formalism, symmetry projectors, tolerances, and disable switches.
- `L`: Liouvillian or Hamiltonian matrix appropriate to the formalism.
- `rho`: initial state for source-state screening, or destination state for destination-state screening. Wavefunctions and Liouville states may form a horizontal stack.

## Reduction path

Input validation rejects nonzero cross-substance blocks with `Spinach:reduce:crossSubstanceGenerator` when `L` has the compiled spin dimension, there is more than one substance, and no reaction records are declared. This check precedes projector construction and applies also to adjoint generators supplied by destination screening. Enlarged spatial-spin inputs are not interpreted using spin-only offsets.

After validation, the source checks whether trajectory-level reduction is disabled by `spin_system.sys.disable` containing `'trajlevel'`; if so, it reports the setting and returns the unit projector `1`. Otherwise the available operations depend on `spin_system.bas.formalism` and the disable settings.

For `zeeman-hilb` and `zeeman-wavef`, the code uses supplied permutation-symmetry irreducible-representation projectors when symmetry treatment is enabled and that information is available. Zero-dimensional irreps are dropped; the state contribution is also screened against `spin_system.tols.irrep_drop`. These formalisms use symmetry screening rather than the Liouville-space zero-track and path-tracing stages.

For `zeeman-liouv` and `sphten-liouv`, the code tries symmetry factorisation when available and not disabled, then applies zero-track elimination when enabled and path tracing to identify disconnected subspaces. Disabling symmetry skips that factorisation; the zero-track and path-tracing stages remain part of this formalism's route, but zero-track elimination returns an identity projector unless `'zte'` is present in `sys.enable`. For compiled spherical-tensor generators, irreps containing substance unit directions survive population screening even when their concentrations are zero. Unit support is mapped into each irrep before zero-track elimination, and included in the subsequent path-tracing screening. Original block offsets are never applied to projected coordinates. The detailed reductions therefore depend on the chosen formalism, supplied symmetry data, input state, and configured tolerances.

`projectors` is a cell array of projectors into independently evolving reduced subspaces. For each `P`, use `L_reduced=P'*L*P` and `rho_reduced=P'*rho` for matrices and state vectors, respectively.

## Header notes

The reduction order is symmetry factorisation, zero-track elimination, then disconnected-subspace identification by path tracing. Further details are in [doi:10.1016/j.jmr.2008.08.008](https://doi.org/10.1016/j.jmr.2008.08.008), [doi:10.1063/1.3398146](https://doi.org/10.1063/1.3398146), and [doi:10.1016/j.jmr.2011.03.010](https://doi.org/10.1016/j.jmr.2011.03.010). See also [the function Wiki page](https://spindynamics.org/wiki/index.php?title=reduce.m).
