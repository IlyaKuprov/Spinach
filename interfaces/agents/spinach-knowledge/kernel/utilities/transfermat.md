# kernel/utilities/transfermat.m

- Signature: `T=transfermat(amp_inps,amp_outs)`

## Purpose

Calculates the transfer matrix for linear filters from paired amplifier input and output vectors.

## Physical / mathematical content

The matrix `T` maps the input vectors to the output vectors in the least-squares sense: `amp_outs = T*amp_inps`.

## Numerical / algorithmic content

The routine checks the input dimensions and obtains the least-squares transfer matrix using MATLAB right matrix division (the source describes this as an SVD pseudoinverse).

## Parameters / inputs

- `amp_inps` — matrix whose columns are amplifier input vectors.
- `amp_outs` — matrix whose columns are the corresponding amplifier output vectors.

## Outputs

- `T` — transfer matrix satisfying `amp_outs = T*amp_inps` in the least-squares sense.

- Note: the number of input-output vector pairs should be bigger than the number of elements in those vectors.

## Implementation structure

- Checks consistency of the input matrices.
- Runs the least-squares solve with `T=amp_outs/amp_inps`.
- Source documentation: <https://spindynamics.org/wiki/index.php?title=transfermat.m>
