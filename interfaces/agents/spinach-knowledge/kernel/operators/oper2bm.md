# kernel/operators/oper2bm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/oper2bm.m`
- Signature: `[states,coeffs]=oper2bm(A)`
- Total lines: 71

## Purpose

Bosonic monomial operator expansion of a user-specified square matrix. Syntax: [states,coeffs]=oper2bm(A)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -a square matrix

## Outputs

- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination

## Implementation structure

- Bosonic monomial operator expansion of a user-specified
- square matrix. Syntax:
- [states,coeffs]=oper2bm(A)
- A -a square matrix
- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `boson_mono()`, `hdot()`, `coeffs()`, `transpose()`, `eps()`, `states()`, `ismatrix()`.
