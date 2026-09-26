# kernel/overloads/@ttclass/ttort.m

- Signature: `[tt,lognrm]=ttort(tt,direct)`

## Purpose

Performs TT-orthogonalisation for a tensor train (or for each tensor train in a buffered sum). Syntax: [tt,lognrm]=ttort(tt,direct)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- direct=+1 -gives left-to-right orthogonality,
- direct=-1 -gives right-to-left orthogonality
- tt -tensor train object, possibly with buffered sums

## Outputs

- tt -tensor train object with all terms in the buffe-
- red sum has all of them orthogonalised in the
- direction requested
- lognrm -if this output is present, all buffered trains
- are also normalized, and natural logs of their
- norms returned in the vector lognrm. Use this
- option if the tensor norm is likely to exceed
- realmax()=1.7977e+308.
- Note: normally you should not call this subroutine directly.

## Implementation structure

- Performs TT-orthogonalisation for a tensor train (or for each tensor
- train in a buffered sum). Syntax:
- [tt,lognrm]=ttort(tt,direct)
- direct=+1 -gives left-to-right orthogonality,
- direct=-1 -gives right-to-left orthogonality
- tt -tensor train object, possibly with buffered sums
- tt -tensor train object with all terms in the buffe-
- red sum has all of them orthogonalised in the
- direction requested
- lognrm -if this output is present, all buffered trains
- are also normalized, and natural logs of their
- norms returned in the vector lognrm. Use this
