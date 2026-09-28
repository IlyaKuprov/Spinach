# kernel/operators/stevens.m

- Signature: `S=stevens(mult,k,q)`

## Purpose

Construct an extended Stevens operator matrix for a spin multiplicity, rank, and projection.

## Physical / mathematical content

- Starts from the spin raising operator raised to rank `k`, then commutes with the lowering operator to obtain the requested projection.
- Uses the Hermitian sum for `q >= 0` and the Hermitian difference divided by `2i` for `q < 0`.

## Numerical / algorithmic content

- Normalization uses explicitly stockpiled integer coefficients. The historical definition is irregular; only ranks up to 12 are available.

## Parameters / inputs

- `mult` — multiplicity of the spin in question.
- `k` — Stevens operator rank, an integer from 0 to 12.
- `q` — Stevens operator projection, an integer from `-k` to `k`.

## Outputs

- `S` — Stevens operator matrix.

## Implementation structure

- Validates the inputs, obtains spin matrices with `pauli(mult)`, constructs and commutes the operator, applies its normalization coefficient, and forms the result according to the sign of `q`.

## Reference

- <https://spindynamics.org/wiki/index.php?title=stevens.m>