# kernel/operators/hilb2liouv.m

- Signature: `L=hilb2liouv(H,conv_type)`
- DIRECT source: [kernel/operators/hilb2liouv.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/hilb2liouv.m)
- Wiki: [hilb2liouv.m](https://spindynamics.org/wiki/index.php?title=hilb2liouv.m)

## Definition and ordering

For a square `d`-by-`d` matrix `H`, the function builds a Liouville-space matrix from Kronecker products with the sparse identity `I=speye(size(H))`. Matrix operators are represented as column-stacked vectors, consistent with MATLAB `H(:)` ordering. In this convention `kron(I,H)` acts by left multiplication, `X -> H*X`, while `kron(transpose(H),I)` acts by right multiplication, `X -> X*H`. The source uses the non-conjugating `transpose(H)`.

## Conversion choices

- `'left'` returns `kron(I,H)`: left action, `X -> H*X`.
- `'right'` returns `kron(transpose(H),I)`: right action, `X -> X*H`.
- `'comm'` returns `kron(I,H)-kron(transpose(H),I)`: `X -> H*X-X*H`.
- `'acomm'` returns `kron(I,H)+kron(transpose(H),I)`: `X -> H*X+X*H`.
- `'statevec'` returns `H(:)`, the column-stacked vector itself, rather than a `d^2`-by-`d^2` superoperator.

For the four action choices and square `H`, `L` is `d^2`-by-`d^2`. The commutator and anticommutator branches contain only the stated difference or sum: there is no extra scalar such as `1i` or `hbar`, and no matrix exponential. They are direct superoperators, not time propagators.

## Inputs and checks

An explicit cell array of per-substance matrices is converted block-wise: action superoperators are assembled with `blkdiag`, while state vectors are stacked vertically. Numeric input retains its original single-matrix meaning; block boundaries are never inferred from zero entries.

The source checks that `H` is numeric or a cell array and `conv_type` is a character array; unrecognised conversion types raise an error. It does not explicitly check that `H` is square, although the matrix-action formulas above presume a square operator.
