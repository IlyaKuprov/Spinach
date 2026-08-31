# kernel/overloads/@ttclass/unit_like.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/unit_like.m`
- Signature: `A=unit_like(A)`
- Total lines: 65

## Purpose

Returns a unit object of the same type as whatever is supplied.

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Unit tensor train of the same topology; implemented by `mode_sizes=sizes(A)`.
- Lines 37-38: Unit sparse matrix of the same dimension; implemented by `A=speye(size(A))`.
- Lines 42-43: Unit dense matrix of the same dimension; implemented by `A=eye(size(A))`.
- Lines 47-48: Complain and bomb out; implemented by `error('the input is not a square matrix or a representation thereof.')`.

### Control flow inferred from the code

- Line 21: conditional branch on `isa(A,'ttclass')`.
- Line 25: conditional branch on `all(mode_sizes(:,1)==mode_sizes(:,2))`.
- Line 27: `for` loop over `k=1:A.ncores`.

### Key state/data transformations

- Lines 24: computes `mode_sizes` using `mode_sizes=sizes(A)`.
- Lines 26: computes `core` using `core=cell(A.ncores,1)`.
- Lines 28: computes `core{k}` using `core{k}=eye(mode_sizes(k,1))`.
- Lines 30: computes `A` using `A=ttclass(1,core,0)`.

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
