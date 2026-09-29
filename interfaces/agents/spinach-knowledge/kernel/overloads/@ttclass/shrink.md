# kernel/overloads/@ttclass/shrink.m

- Signature: `ttrain=shrink(ttrain)`

## Action

`shrink` first calls `pack` to absorb the buffered sum of trains into one TT, then calls `ttort(ttrain,+1)` for left-to-right orthogonalisation. It forms `nrm=ttrain.coeff*norm(ttrain.cores{d,1}(:),2)`; if this value is exactly zero, it returns `0*unit_like(ttrain)` without truncating. Otherwise it calls `truncate` to reduce TT ranks according to the train's truncation settings. The core count and physical mode dimensions are retained; ranks may be reduced. When more than one train is packed, `pack` incorporates each train's coefficient into its first core and sets the packed coefficient to one; its single-train fast path leaves the input coefficient as-is. Subsequent coefficient handling is delegated to orthogonalisation and truncation.

The result remains a TT, rather than a materialised matrix, except when every core has singleton physical dimensions in both axes (core dimensions 2 and 3); in that case the method calls `full` and returns the scalar value. This test concerns physical dimensions, not bond dimensions. The method has no explicit conjugation or transpose action.

## Input and output

- `ttrain` — tensor-train object.
- `ttrain` — packed, left-to-right orthogonalised, truncated train, or the scalar result in the singleton-physical-mode case.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/shrink.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/shrink.m)

D. Savostyanov and I. Kuprov.
