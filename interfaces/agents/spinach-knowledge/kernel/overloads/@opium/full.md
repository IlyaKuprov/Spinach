# kernel/overloads/@opium/full.m

- Signature: `M=full(M)`

## Meaning and behaviour

An `opium` object represents the scaled unit matrix `coeff*I_dim`. This overload returns `M.coeff*eye(M.dim)`: a square numeric matrix, with the identity explicitly formed by `eye`. It therefore materialises the represented matrix rather than retaining the compact `opium` object. This source contains no explicit input validation.

## Input and output

- Input: `M`, an `opium` object with the `coeff` and `dim` fields used by the implementation.
- Output: the full scaled unit matrix of order `M.dim`.

## Source links

- MATLAB source: [kernel/overloads/@opium/full.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/full.m)
- Existing Wiki page: [opium/full.m](https://spindynamics.org/wiki/index.php?title=opium/full.m)
