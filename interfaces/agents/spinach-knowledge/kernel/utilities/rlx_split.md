# kernel/utilities/rlx_split.m

- Signature: `[R1,R2,Rm]=rlx_split(spin_system,R)`

## Purpose

Partitions a relaxation superoperator into longitudinal single-spin, transverse single-spin, and residual components. The routine is available only in the `sphten-liouv` formalism.

## Physical / mathematical content

The longitudinal block acts within basis states representing purely longitudinal single-spin orders; the transverse block acts within purely transverse single-spin orders. All matrix elements not included in those two blocks are retained in the residual component.

## Numerical / algorithmic content

The function maps the basis labels to angular-momentum ranks `L` and projections `M`. Rows with exactly one nonzero spin rank are single-spin orders. A state is classified as longitudinal when its nonzero rank has `M=0`, or transverse when it has `M~=0`. The two diagonal submatrices of `R` are extracted and `Rm=R-R1-R2`.

## Parameters / inputs

- `spin_system` - Spinach system structure whose basis formalism is `sphten-liouv`.
- `R` - square numeric relaxation superoperator in that basis.

## Outputs

- `R1` - part of `R` acting within purely longitudinal single-spin states.
- `R2` - part of `R` acting within purely transverse single-spin states.
- `Rm` - remainder `R-R1-R2`, including mixed and other matrix elements.

## Implementation structure

The routine validates the basis formalism and that `R` is a square matrix, computes the basis rank/projection labels with `lin2lm`, constructs the longitudinal and transverse single-spin masks, and places each selected block into a zero matrix of the same size as `R`.

## Reference

[Spin Dynamics Wiki: rlx_split.m](https://spindynamics.org/wiki/index.php?title=rlx_split.m)
