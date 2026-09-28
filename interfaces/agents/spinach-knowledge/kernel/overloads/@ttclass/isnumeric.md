# kernel/overloads/@ttclass/isnumeric.m

- Signature: `answer=isnumeric(tt)`

## Purpose

Reports whether the input is a non-empty tensor-train object according to the `ttclass` representation.

## Numerical / algorithmic content

Returns true only when `tt` is a `ttclass` object and `tt.cores` is non-empty. Other inputs return false; this predicate checks the core container, not the numerical values or whether the represented tensor is nonzero.

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true for non-empty tensor train objects

## Implementation structure

The implementation combines `isa(tt,'ttclass')` with `~isempty(tt.cores)` and assigns the resulting logical flag.
