# kernel/grids/grid_test.m

- Signature: `grid_profile=grid_test(alphas,betas,gammas,weights,ranks,sfun)`

## Behaviour

For each requested spherical rank `l`, forms the weighted Wigner matrix `D=sum(weights(j)*wigner(l,alphas(j),betas(j),gammas(j)))`. Each matrix is `(2*l+1)` by `(2*l+1)`; the returned `grid_profile` has the same shape as `ranks`. Euler-angle inputs are finite real column vectors in radians, with one entry per grid point; `alphas` may be zero for single-angle grids. `weights` is a matching column vector of finite positive real values. The code uses weights as supplied and does not check that they sum to one or renormalise them.

`ranks` is a vector of finite nonnegative integers. The selector `sfun` chooses the reported statistic: `'D_lmn'` subtracts `krondelta(0,l)` from the spectral 2-norm of the whole matrix; `'Y_lm'` subtracts it from the 2-norm of the central row `D(l+1,:)`; `'Y_l0'` subtracts it from the central element `D(l+1,l+1)`. These are the three-angle, two-angle, and single-angle diagnostics respectively. The returned values are the code's norm-minus-delta scores, not absolute-valued errors. Each score is also reported; if no output is requested, the function plots the profile against rank.

Angles are in radians; the source assigns no physical unit to `weights`. No time, frequency-offset, eigenfield, or evolution input is part of this grid diagnostic.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_test.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grid_test.m)