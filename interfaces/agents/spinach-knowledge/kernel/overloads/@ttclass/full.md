# kernel/overloads/@ttclass/full.m

- Signature: `answer=full(ttrain)`

## Purpose

Converts a tensor-train representation of a matrix into a dense matrix.

## Input

- `ttrain` — tensor train object.

## Output

- `answer` — full matrix represented by `ttrain`.

## Note

The dense result can be very large; careless use may exhaust available memory.

## Algorithm

The function preallocates the result, obtains the train ranks and mode sizes, contracts the cores from the last core toward the first for each train, and adds each contracted matrix weighted by its coefficient.
