# kernel/kinetics/equilibrate.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/kinetics/equilibrate.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=equilibrate.m)

`equilibrate(K,c0)` returns the steady-state concentration vector for the linear chemical-kinetics equation `dc/dt = K*c`. It is not a line-shape routine: there is no frequency axis, so Hz-versus-angular-frequency conversion does not apply.

The solution is constrained by `K*c = 0` and `sum(c) = sum(c0)`. The code solves the stacked system `[ones(1,n); K]*c = [sum(c0); zeros(n,1)]`, where `n` is the number of species. It partitions disconnected reaction components using the nonzero pattern of `K` or its transpose and solves each component recursively. Before solving, each nonzero component's `K` is divided by `max(abs(K(:)))`. A common rescaling of all rates leaves the mathematical equilibrium unchanged, provided the absolute column-sum guard still passes; use one consistent inverse-time unit (for example, s^-1). There is no frequency axis, and the absolute rate scale is not returned.

## Inputs and outputs

`K` must be a real numeric square matrix whose column sums each have absolute value strictly less than `10*eps('double')`. This is the implemented mass-conservation guard. `c0` must be a real numeric non-negative column vector with one entry per column of `K`. The routine does not separately enforce kinetic-generator sign conditions on `K`'s off-diagonal entries or check that the solved concentrations are non-negative.

The output `c` is an `n`-by-1 concentration vector with the same total as `c0`. If `c0` is the zero vector, it returns `c0` directly. Otherwise it errors if the stacked system's condition number exceeds `1/sqrt(eps('double'))`; below that threshold it solves with MATLAB backslash.
