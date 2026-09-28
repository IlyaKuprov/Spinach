# kernel/overloads/@ttclass/ismatrix.m

- Signature: `answer=ismatrix(tt)`

## Purpose

Tests whether the input is a non-empty tensor-train object.

## Input

- `tt` — object to test.

## Output

- `answer` — logical true when `tt` is a `ttclass` object with non-empty cores; logical false otherwise.

## Behavior

The function checks the object's class and whether its `cores` property is non-empty.
