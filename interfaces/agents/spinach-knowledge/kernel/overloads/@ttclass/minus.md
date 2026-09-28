# kernel/overloads/@ttclass/minus.m

- Signature: `a=minus(a,b)`

## Purpose

Represents the difference of two tensor trains without immediately performing a recompression.

## Numerical / algorithmic content

The inputs must be `ttclass` objects with matching core counts and mode sizes. The result is assembled by concatenating their train components and pairing the second operand's coefficients with negative signs. Zero-coefficient components are discarded; if none remain, the result is set to a zero tensor train.

## Parameters / inputs

- a -a tensor train object
- b -a tensor train object

## Outputs

- a -the tensor train representing the difference `a-b`
