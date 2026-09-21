# kernel/overloads/@ttclass/kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/kron.m`
- Signature: `c=kron(a,b)`
- Total lines: 67

## Purpose

Kronecker product of two matrices in a tensor train format. Syntax: c=kron(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a,b -tensor train objects

## Outputs

- c -a tensor trin object
- WARNING: the result is not the same as the flat matrix Kronecker pro-
- duct (it is a row and column permutation away from it), but
- the resulting order of elements is consistent with the out-
- put of the tensor train vectorization (ttclass/vec) operati-
- on output.

## Implementation structure

- Kronecker product of two matrices in a tensor train format. Syntax:
- c=kron(a,b)
- a,b -tensor train objects
- c -a tensor trin object
- WARNING: the result is not the same as the flat matrix Kronecker pro-
- duct (it is a row and column permutation away from it), but
- the resulting order of elements is consistent with the out-
- put of the tensor train vectorization (ttclass/vec) operati-
- on output.
- Shrink a and b before going any further
- Read sizes and ranks of the operands
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `shrink()`, `ranks()`, `sizes()`, `a_ranks()`, `a_sizes()`, `b_ranks()`, `b_sizes()`.
