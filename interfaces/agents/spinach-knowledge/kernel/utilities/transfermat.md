# kernel/utilities/transfermat.m

## Purpose

Computes the transfer matrix of a linear filter from stacks of observed amplifier input and output vectors ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/transfermat.m)).

## Behaviour

- Syntax: `T=transfermat(amp_inps,amp_outs)`.
- The function first runs a consistency check (`grumble`) on the two input stacks, then computes `T=amp_outs/amp_inps`, described in the header as the SVD pseudoinverse route.
- The returned matrix satisfies `amp_outs=T*amp_inps` in the least squares sense.
- The header notes that the number of input-output vector pairs should be bigger than the number of elements in those vectors.
- Consistency enforcement (`grumble`) errors when:
  - `amp_inps` is not numeric, or has fewer columns than rows (`size(amp_inps,2)<size(amp_inps,1)`), with message `amp_inps must be a stack of column vectors wider than it is tall.`;
  - `amp_outs` is not numeric, or has fewer columns than rows, with message `amp_outs must be a stack of column vectors wider than it is tall.`;
  - the two stacks have different numbers of vectors (`size(amp_inps,2)~=size(amp_outs,2)`), with message `the number of vectors in amp_inps and amp_outs stacks must be the same.`.

## Inputs and outputs

Inputs:

- `amp_inps` — numeric matrix with amplifier input vectors as columns; must have at least as many columns as rows.
- `amp_outs` — numeric matrix with amplifier output vectors as columns; must have at least as many columns as rows and the same number of columns as `amp_inps`.

Outputs:

- `T` — the transfer matrix, such that `amp_outs=T*amp_inps` in the least squares sense.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/transfermat.m>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=transfermat.m>
