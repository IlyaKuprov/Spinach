# kernel/overloads/@ttclass/unit_like.m

- Signature: `A=unit_like(A)`

## Purpose

Returns a unit object of the same type as whatever is supplied.

## Physical / mathematical content

Returns an identity matrix or tensor-train representation for a square matrix.

## Numerical / algorithmic content

For a tensor-train input, the function checks that each core has matching row and column mode sizes, creates an identity matrix for each core, and constructs `ttclass(1,core,0)`. For a sparse square matrix it returns `speye(size(A))`; for a dense square matrix it returns `eye(size(A))`. Other inputs, including tensor trains that do not represent square matrices, raise an error.

## Syntax

```matlab
A=unit_like(A)
```

## Parameters / inputs

- `A` — a full or sparse square matrix, or a tensor-train representation of a square matrix.

## Outputs

- `A` — a unit matrix in the same format.
