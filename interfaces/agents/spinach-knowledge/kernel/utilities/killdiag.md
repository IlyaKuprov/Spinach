# kernel/utilities/killdiag.m

## Purpose

`killdiag.m` zeroes out a band along the diagonal of a 2D spectrum using a brush of a specified width. It is documented as part of the Spinach kernel utilities ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/killdiag.m)).

## Behaviour

- Syntax: `spec=killdiag(spec,brush_dim)`.
- For each column `n` of the spectrum, the row index on the diagonal is computed as `k=n*size(spec,1)/size(spec,2)`.
- The brush extends from `round(k-(brush_dim-1)/2)` to `round(k+(brush_dim-1)/2)`, giving a band of `brush_dim` points centred on the diagonal.
- Row indices outside the array boundaries (`k<1` or `k>size(spec,1)`) are discarded before zeroing.
- The selected elements `spec(k,n)` are set to zero; the function returns the modified matrix.
- Input consistency is enforced by the internal `grumble` function, which errors when: `spec` is not a numeric matrix; `brush_dim` is not a positive real integer scalar; or the brush is wider than the spectrum (`any(size(spec)<brush_dim)`).

## Inputs and outputs

**Inputs**

- `spec` — 2D matrix representing a spectrum.
- `brush_dim` — the width of the band to zero out around the diagonal, in points.

**Outputs**

- `spec` — 2D matrix representing the spectrum with the diagonal band zeroed.

## References

- Spin Dynamics Wiki page for `killdiag.m`: <https://spindynamics.org/wiki/index.php?title=killdiag.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/killdiag.m>
