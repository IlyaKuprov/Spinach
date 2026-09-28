# kernel/overloads/@ttclass/dot.m

- Signature: `c=dot(a,b)`

## Purpose

Computes the inner product of two tensor-train representations of numerical arrays.

## Inputs

- `a`, `b` — tensor-train objects with matching mode sizes and the same number of cores.

## Output

- `c` — inner product of `a` and `b`.

## Algorithm

The function checks that both inputs are tensor trains with compatible sizes, then evaluates the product as `ctranspose(a)*b`.
