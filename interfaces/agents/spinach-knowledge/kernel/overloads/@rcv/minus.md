# kernel/overloads/@rcv/minus.m

- Signature: `A=minus(A,B)`

## Purpose

Computes `A-B` for two RCV operands by negating `B` and delegating the addition to `plus`.

## Mathematical content

The operation is subtraction of the right operand from the left operand; this function implements it as `A=plus(A,(-1)*B)`.

## Numerical / algorithmic content

The function calls its consistency check before dispatching to scalar multiplication and `plus`. The check raises an error only when neither operand is an RCV object.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- A -result A-B as an RCV sparse matrix

## Implementation structure

- Run the local input check.
- Replace `A` with `plus(A,(-1)*B)` and return the result.
