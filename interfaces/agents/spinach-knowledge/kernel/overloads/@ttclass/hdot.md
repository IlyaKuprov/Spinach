# kernel/overloads/@ttclass/hdot.m

- Signature: `c=hdot(a,b)`

## Purpose

Computes the Hadamard dot product of two tensor-train matrices.

## Inputs

- `a`, `b` — tensor-train objects representing arrays with the same dimensions and internal topology.

## Output

- `c` — scalar Hadamard dot product of `a` and `b`.

## Behavior

The function requires both inputs to be `ttclass` objects with the same number of cores and matching mode sizes. For each pair of train-buffer entries, it reshapes corresponding core data, forms `core_a' * core_b` while iterating over cores, and adds the resulting `x` to `c`.
