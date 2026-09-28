# kernel/overloads/@ttclass/conj.m

- Signature: `tt=conj(tt)`

## Purpose

Complex-conjugates the tensor train by conjugating every core element and every coefficient.

## Input

- `tt` — tensor train object.

## Output

- `tt` — tensor train object with conjugated cores and coefficients.

## Algorithm

The function loops over all trains and cores, applies `conj` to each core, then applies `conj` to the coefficient array.
