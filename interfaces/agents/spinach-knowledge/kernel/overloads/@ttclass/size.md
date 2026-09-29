# kernel/overloads/@ttclass/size.m

- Signature: `varargout=size(tt,dim)`

## Action

The method multiplies each core's second dimension across the core sequence to obtain the represented row count `m`, and each third dimension to obtain the column count `n`; it uses the first train's products. Supported forms are `sz=size(tt)` (returns `[m n]`), `[m,n]=size(tt)`, and `d=size(tt,dim)` for `dim=1` or `dim=2`. A supplied dimension other than 1 or 2, or an unsupported input/output form, reaches the `incorrect call syntax.` error. If either returned dimension exceeds MATLAB's `intmax`, it raises the tensor-train-dimensions error.

This is a metadata query: it does not change cores or ranks and does not materialise the represented matrix. It applies no conjugation or transpose.

## Input and output

- `tt` — tensor-train representation of a matrix.
- `dim` — optional selector: 1 for rows or 2 for columns.
- Output — row and column dimensions, separately or as a two-element row vector.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/size.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/size.m)

D. Savostyanov and I. Kuprov.
