# kernel/overloads/@ttclass/unit_like.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/unit_like.m`
- Signature: `A=unit_like(A)`
- Total lines: 65

## Purpose

Returns a unit object of the same type as whatever is supplied.

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Syntax

```matlab
A=unit_like(A)
```

## Parameters / inputs

- A -a full or sparse square matrix, or a tensor train
- representation of a square matrix

## Outputs

- A -a unit matrix in the same format

## Implementation structure

- Returns a unit object of the same type as whatever is supplied.
- A=unit_like(A)
- A -a full or sparse square matrix, or a tensor train
- representation of a square matrix
- A -a unit matrix in the same format
- Unit tensor train of the same topology
- Unit sparse matrix of the same dimension
- Unit dense matrix of the same dimension
- Complain and bomb out
- Briefly stated, the Gell-Mann Amnesia effect is as follows. You open the
- newspaper to an article on some subject you know well. You read the arti-
- cle and see the journalist has absolutely no understanding of either the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sizes()`, `all()`, `mode_sizes()`, `ttclass()`, `ismatrix()`, `issparse()`, `speye()`.
