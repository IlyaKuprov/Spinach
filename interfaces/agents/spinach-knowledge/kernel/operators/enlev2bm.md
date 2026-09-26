# kernel/operators/enlev2bm.m

- Signature: `[states,coeffs]=enlev2bm(nlevels,lvl_num)`

## Purpose

Bosonic monomial expansion of specific bosonic energy level projectors. Syntax: [states,coeffs]=enlev2bm(nlevels,lvl_num)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- nlevels -number of energy levels in the mode,
- a positive integer
- lvl_num -energy level number, counting from the
- empty mode state upwards

## Outputs

- states -states, in the Spinach BM basis index-
- ing, that contribute to the operator in
- question; use lin2kq to convert to K,Q
- bosonic monomial indices
- coeffs -coefficients with which the BMs enter
- the linear combination

## Implementation structure

- Bosonic monomial expansion of specific bosonic energy
- level projectors. Syntax:
- [states,coeffs]=enlev2bm(nlevels,lvl_num)
- nlevels -number of energy levels in the mode,
- a positive integer
- lvl_num -energy level number, counting from the
- empty mode state upwards
- states -states, in the Spinach BM basis index-
- ing, that contribute to the operator in
- question; use lin2kq to convert to K,Q
- bosonic monomial indices
- coeffs -coefficients with which the BMs enter
