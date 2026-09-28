# kernel/overloads/@ttclass/sum.m

- Signature: `answer=sum(ttrain,dim)`

## Purpose

Sum the elements of a tensor-train matrix along its row or column dimension.

## Parameters / inputs

- `ttrain` — tensor train representing a matrix.
- `dim` — summation dimension, 1 or 2; when omitted, the first non-singleton dimension is selected, matching MATLAB's matrix behavior.

## Outputs

- `answer` — tensor train representing the summation result, or its full scalar value when all resulting modes are singleton.

## Implementation

If all physical dimensions are singleton, the function immediately returns `full(ttrain)`. Otherwise, it creates an auxiliary train with the same coefficients and zero tolerances, sums each core over the selected physical dimension, and reshapes that mode to size one while retaining the other physical dimension and the TT ranks. If the resulting train has only singleton modes, it is converted to a scalar with `full`.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/sum.m`](https://spindynamics.org/wiki/index.php?title=ttclass/sum.m).
