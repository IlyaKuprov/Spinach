# kernel/conventions/transforms/tsm2param.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/tsm2param.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=tsm2param.m)

## Purpose and limitation

Converts a traceless symmetric interaction matrix to axiality, rhombicity, and Euler angles in Mehring eigenvalue order. The source warns that the parameterisation is unstable and recommends publishing the 3x3 matrix instead, consistent with IUPAC guidance.

## Input and reconstruction

`M` is real numeric input with either nine matrix elements or five independent elements ordered `[Mxx, Mxy, Mxz, Myy, Myz]`. For five elements the function constructs:

`[Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz -Mxx-Myy]`

For nine elements, the source checks `issymmetric(M)` and `abs(trace(M))<=10*eps`; it has no explicit `size(M)==[3 3]` or finite-value check. Inputs passing those checks are passed to `eig`, which requires a square matrix. The five-element input is indexed linearly in the stated order.

## Conversion and outputs

Let `lambdaX` be the smallest eigenvalue, `lambdaZ` the largest, and `lambdaY` the remaining eigenvalue. The returned scalar invariants are `ax=2*lambdaZ-(lambdaX+lambdaY)` and `rh=lambdaY-lambdaX`. They retain the input matrix's units; no unit conversion is applied.

The eigenvector columns are reordered as X, Y, Z and multiplied by `det(V)` so the resulting orientation matrix has determinant +1. The function passes that matrix to [`dcm2euler.m`](dcm2euler.md); `angles` is its 1-by-3 ZYZ active Euler-angle row in radians. The eigensystem does not define a unique orientation when eigenvalues are degenerate, consistent with the source's warning about instability.

The implementation calls its validation helper before diagonalisation. It requires real numeric input with five or nine elements; on the nine-element path it additionally checks symmetry and the absolute trace tolerance above.
