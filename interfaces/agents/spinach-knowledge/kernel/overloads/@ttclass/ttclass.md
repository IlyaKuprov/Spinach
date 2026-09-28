# kernel/overloads/@ttclass/ttclass.m

- Signature: `tt=ttclass(coeff,kronterms,tolerance)`

## Purpose

Creates an object of a tensor train class. A tensor train is a type of un-opened Kronecker product that behaves as a matrix or a vector of a very large dimension, but takes a reasonable amount of memory to store. See https://doi.org/10.1137/090752286 for further information. Syntax: tt=ttclass(coeff,kronterms,tolerance)

## Physical / mathematical content

A tensor train represents a matrix or vector assembled from Kronecker-product terms in a compact form. The cited reference is https://doi.org/10.1137/090752286.

## Numerical / algorithmic content

The constructor takes coefficients, Kronecker terms, and tolerances. When multiple columns are supplied, the header notes describe the result as the sum of the corresponding individual tensor trains.

## Parameters / inputs

- `coeff` — coefficient in front of the spin operator, usually the interaction magnitude.
- `kronterms` — column cell array of matrices whose Kronecker product makes up the spin operator.
- `tolerance` — maximum deviation in the 2-norm between the tensor-train representation and the flat matrix representation that the tensor-train format is allowed to introduce.

## Outputs

- `tt` — tensor-train object.

## Header notes

1. If multiple columns are supplied in `kronterms`, multiple coefficients are given in `coeff`, and multiple tolerances are given in `tolerance`, the resulting tensor train is assumed to be the sum of the individual tensor trains specified in different columns.
2. Tensor trains are exotic and capricious structures; do not use them unless you know what you are doing.
