# kernel/reduce.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/reduce.m`
- Signature: `projectors=reduce(spin_system,L,rho)`
- Total lines: 311

## Purpose

Symmetry and trajectory-level state-space reduction. Tries applicable reduction methods unless disabled and returns projectors into independently evolving reduced subspaces.

## Physical / mathematical content

Permutation symmetry separates irreducible sectors; zero-track elimination removes coordinates unoccupied by the trajectory; path tracing separates disconnected subspaces. The available reductions depend on the formalism and the system's disable flags.

## Numerical / algorithmic content

Horizontal stacks of wavefunctions or Liouville states retain their actual complex columns during symmetry screening. The matrix 1-norm of a projected stack is its largest column 1-norm, so occupancy is tested without averaging columns or discarding their phases. Liouville stacks continue through ZTE and path tracing without disabling useful reduction. ZTE streams columns independently and combines their row maxima, avoiding dense whole-stack screening storage and cross-column scaling.

Hilbert-space cell arrays of density matrices retain their existing absolute-value representative for screening. Input validation precedes reduction.

## Syntax

```matlab
projectors=reduce(spin_system,L,rho)
```

## Parameters / inputs

- `spin_system`: Spinach system, including formalism, symmetry projectors, tolerances, and disable switches.
- `L`: Liouvillian or Hamiltonian matrix appropriate to the formalism.
- `rho`: initial state for source-state screening, or destination state for destination-state screening. Wavefunctions and Liouville states may form a horizontal stack.

## Outputs

`projectors` is a cell array of projectors into independently evolving reduced subspaces. For each `P`, use `L_reduced=P'*L*P` and `rho_reduced=P'*rho` for matrices and state vectors, respectively.

## Header notes

The reduction order is symmetry factorisation, zero-track elimination, then disconnected-subspace identification by path tracing. Further details are in [doi:10.1016/j.jmr.2008.08.008](https://doi.org/10.1016/j.jmr.2008.08.008), [doi:10.1063/1.3398146](https://doi.org/10.1063/1.3398146), and [doi:10.1016/j.jmr.2011.03.010](https://doi.org/10.1016/j.jmr.2011.03.010). See also [the function Wiki page](https://spindynamics.org/wiki/index.php?title=reduce.m).
