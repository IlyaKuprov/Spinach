# kernel/overloads/@ttclass/clearcoeff.m

- Signature: `tt=clearcoeff(tt)`

## Purpose

Distributes each tensor-train coefficient across its cores and sets the coefficient to one, without changing the represented tensor train.

## Input

- `tt` — tensor train object.

## Output

- `tt` — tensor train object with each coefficient distributed into its cores and the coefficient array set to one.

## Algorithm

For each train, the function takes the `ncores`-th root of its coefficient, multiplies every core in that train by this factor, then sets that train's coefficient to one.
