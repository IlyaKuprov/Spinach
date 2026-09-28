# kernel/operators/twospinist.m

- Signature: `T=twospinist(spin_system,spin_a,spin_b,indices,type)`

## Purpose

Construct a two-spin irreducible spherical tensor operator.

## Physical / mathematical content

- `indices=[L,M]` selects rank `L=1` or `L=2` and its projection `M`. The implementation provides `M=-1,0,+1` for rank 1 and `M=-2,-1,0,+1,+2` for rank 2.
- Each tensor is assembled from products of raising, lowering, and longitudinal spin operators on the two specified spins.

## Numerical / algorithmic content

- The function selects the requested tensor expression by rank and projection, constructs its terms with `operator(...,type,'csc')`, and combines them with the specified coefficients.
- Invalid rank or projection selections produce an error.

## Parameters / inputs

- `spin_system` — spin system passed to `operator`.
- `spin_a` — number of the first spin.
- `spin_b` — number of the second spin; it must differ from `spin_a`. Both spin numbers must be integers from 1 through the number of spins in the system.
- `indices` — two-element real vector `[L,M]`, with `L=1` or `L=2` and integer `M`.
- `type` — character string. In Liouville space, `'left'` produces a left-side product superoperator, `'right'` a right-side product superoperator, `'comm'` a commutation superoperator (default), and `'acomm'` an anticommutation superoperator. In Hilbert space, `type` is ignored and the operator itself is returned.

## Outputs

- `T` — irreducible spherical tensor operator.

## Implementation structure

- A consistency check validates the spin numbers, `type`, and tensor indices. Nested rank and projection switches then construct the selected tensor; unsupported selections raise an error.

<https://spindynamics.org/wiki/index.php?title=twospinist.m>