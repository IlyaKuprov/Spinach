# kernel/overloads/@ttclass/shrink.m

- Signature: `ttrain=shrink(ttrain)`

## Purpose

Approximate a tensor train with lower TT ranks by orthogonalising it and applying the train's truncation settings.

## Parameters / inputs

- `ttrain` — tensor train object.

## Outputs

- `ttrain` — compressed tensor train, left-to-right orthogonalised; converted to a scalar representation when all core bond dimensions are one.

## Implementation

The function packs the train, applies left-to-right orthogonalisation with `ttort(ttrain,+1)`, and checks its norm. If the norm is zero, it returns a zero-like train immediately. Otherwise it calls `truncate`; if every core has singleton left and right bond dimensions, it converts the result with `full`.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/shrink.m`](https://spindynamics.org/wiki/index.php?title=ttclass/shrink.m).
