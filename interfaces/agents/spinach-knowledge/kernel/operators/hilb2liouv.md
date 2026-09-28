# kernel/operators/hilb2liouv.m

- Signature: `L=hilb2liouv(H,conv_type)`

## Purpose

Converts a Hilbert-space operator into a Liouville-space superoperator or, for `statevec`, a column-stacked state vector.

## Physical / mathematical content

The conversion type selects left multiplication, right multiplication, a commutator, an anticommutator, or direct vectorization of `H`.

## Numerical / algorithmic content

With `I=speye(size(H))`, the returned matrices are formed as follows:

- `left`: `kron(I,H)`.
- `right`: `kron(transpose(H),I)`.
- `comm`: `kron(I,H)-kron(transpose(H),I)`.
- `acomm`: `kron(I,H)+kron(transpose(H),I)`.
- `statevec`: `H(:)`, using MATLAB column-major ordering.

## Parameters / inputs

- H - numeric Hilbert-space operator.
- conv_type - character string selecting the conversion: `'left'`, `'right'`, `'comm'`, `'acomm'`, or `'statevec'`.

## Outputs

- L - resulting Liouville-space superoperator or column-stacked state vector.

## Implementation structure

1. Check that `H` is numeric and `conv_type` is a character string.
2. Create the sparse identity with the dimensions of `H`.
3. Apply the Kronecker-product formula for the selected conversion, or return `H(:)` for `statevec`. Unknown conversion strings raise an error.
