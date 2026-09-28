# kernel/operators/mprealloc.m

- Signature: `A=mprealloc(spin_system,nnzpc)`

## Purpose

Preallocates an all-zero sparse operator of the dimension required by the current Spinach formalism, with storage estimated from the expected number of nonzeros per column.

## Physical / mathematical content

The matrix dimension is taken from the current basis for `sphten-liouv`, from the product of spin multiplicities for `zeeman-wavef` and `zeeman-hilb`, and from the product of squared multiplicities for `zeeman-liouv`.

## Numerical / algorithmic content

The routine calls `spalloc(problem_dim,problem_dim,nnzpc*problem_dim)` to create the sparse matrix. It uses `spin_system.bas.basis` to obtain the `sphten-liouv` dimension and `spin_system.comp.mults` for the Zeeman formalisms. Other formalism values raise an error.

## Parameters / inputs

- spin_system - Spinach system structure containing the current formalism; it must include `bas.formalism`, and the selected case also uses the basis or multiplicities described above.
- nnzpc - expected number of nonzeros per column; the interface describes this as a positive real integer.

## Outputs

- A - all-zero sparse matrix sized for the current formalism and preallocated for `nnzpc*problem_dim` nonzero entries.

## Implementation structure

1. Check the formalism field and the numeric, real, scalar, integer form of `nnzpc`.
2. Select the matrix dimension for `sphten-liouv`, `zeeman-wavef`, `zeeman-hilb`, or `zeeman-liouv`.
3. Allocate the square sparse matrix with `spalloc`; unsupported formalism values raise an error.
