# kernel/overloads/@ttclass/isreal.m

- Signature: `answer=isreal(tt)`

## Purpose

Returns TRUE for real-valued tensor train objects. Syntax: answer=isreal(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true when all coefficients and core
- elements of the tensor train are real

## Implementation structure

- Returns TRUE for real-valued tensor train objects. Syntax:
- answer=isreal(tt)
- tt -tensor train object
- answer -logical true when all coefficients and core
- elements of the tensor train are real
- Non-empty tensor trains should return true()
- Check coefficient first
- If the coefficients are real, check the cores
- Complain and bomb out
- Democracy is a pathetic belief in the collective wisdom
- of individual ignorance.
- H.L. Mencken
