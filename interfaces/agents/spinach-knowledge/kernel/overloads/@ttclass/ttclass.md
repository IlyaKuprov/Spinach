# kernel/overloads/@ttclass/ttclass.m

- Signature: `tt=ttclass(coeff,kronterms,tolerance)`

## Purpose

Creates an object of a tensor train class. A tensor train is a type of un-opened Kronecker product that behaves as a matrix or a vector of a very large dimension, but takes a reasonable amount of memory to store. See https://doi.org/10.1137/090752286 for further informa- tion. Syntax: tt=ttclass(coeff,kronterms,tolerance)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- coeff -coefficient in front of the spin operator,
- usually the interaction magnitude
- kronterms -column cell array of matrices whose Krone-
- cker product makes up the spin operator
- tolerance -maximum deviation in the 2-norm between the
- TT representation and the flat matrix repre-
- sentation that the TT format is allowed to
- introduce

## Outputs

- tt -tensor train object

## Header notes

- 1. If multiple columns are supplied in kronterms, multiple coeffi-
- cients are given in coeff, and multiple tolerances are given in
- tolerance, the resulting tensor train is assumed to be the sum
- of the individual tensor trains specified in different columns.
- 2. Tensor trains are exotic and capricious structures, do not use
- them unless you know what you are doing.

## Implementation structure

- Creates an object of a tensor train class. A tensor train is a type
- of un-opened Kronecker product that behaves as a matrix or a vector
- of a very large dimension, but takes a reasonable amount of memory
- to store. See https://doi.org/10.1137/090752286 for further informa-
- tion. Syntax:
- tt=ttclass(coeff,kronterms,tolerance)
- coeff -coefficient in front of the spin operator,
- usually the interaction magnitude
- kronterms -column cell array of matrices whose Krone-
- cker product makes up the spin operator
- tolerance -maximum deviation in the 2-norm between the
- TT representation and the flat matrix repre-
