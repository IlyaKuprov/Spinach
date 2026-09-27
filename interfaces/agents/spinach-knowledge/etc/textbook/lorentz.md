# etc/textbook/lorentz.m

- Signature: `[J,K,Kil]=lorentz(L)`

## Purpose

Constructs matrix generators for the direct-sum Lorentz-group representation `(L,0) ⊕ (0,L)`, with inversion, for positive integer or half-integer rank `L`.

## Construction

Let `D=2L+1` and let `s` be the spin-`L` matrices returned by `pauli(D)`. The rotation generators are block diagonal with `s` in both blocks; the boost generators have `+i s` in the first block and `-i s` in the second. Each returned generator is a `2D × 2D` matrix.

If the third output is requested, the function also computes the 6-by-6 Killing form. It forms the adjoint-representation matrices for the six generators and evaluates their pairwise trace products; this calculation is the expensive optional part.

## Inputs

- `L` — positive integer or half-integer representation rank, supplied as a real numeric scalar.

## Outputs

- `J` — structure containing rotation generators `J.x`, `J.y`, and `J.z`.
- `K` — structure containing boost generators `K.x`, `K.y`, and `K.z`.
- `Kil` — 6-by-6 Killing form; computed only when requested.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=lorentz.m).
