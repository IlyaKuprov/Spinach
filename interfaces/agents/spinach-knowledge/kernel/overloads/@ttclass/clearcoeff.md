# kernel/overloads/@ttclass/clearcoeff.m

- Signature: `tt=clearcoeff(tt)`

## Purpose

Absorbs physical coefficients into tensor train cores without changing the value represented by the tensor train. Syntax: tt=clearcoeff(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object

## Outputs

- tt -tensor train object with each coefficient distributed
- into its cores and the coefficient array set to one

## Implementation structure

- Absorbs physical coefficients into tensor train cores without
- changing the value represented by the tensor train. Syntax:
- tt=clearcoeff(tt)
- tt -tensor train object
- tt -tensor train object with each coefficient distributed
- into its cores and the coefficient array set to one
- Get the number of cores and trains
- Loop over the trains in the buffer
- Scale the coefficient
- Apply it to cores
- Erase the coefficient
- "Moral outrage is a middle-class luxury."
