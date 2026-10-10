# kernel/grids/grid_fibon.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_fibon.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=grid_fibon.m)

## Purpose

Build a deterministic two-angle quadrature grid on the sphere. The cited construction is Appendix A.5 of [the original paper](http://dx.doi.org/10.1016/j.jmr.2014.05.009).

## Inputs and outputs

- `type` must be a character array; `parm` must be a positive real integer.
- For a selected grid with N points, `alps`, `bets`, and `gams` are N-by-1 Euler-angle columns in radians. `alps` is zero throughout because these grids use two angles.
- `whts` contains one normalised spherical-area weight per point; `vorn` contains the corresponding Voronoi tessera vertex data.

## Grid rules

The executable switch accepts `'fib'`, `'zcw'`, and `'zcwn'`. For `'fib'`, with n = `parm`, k runs from -n to n, beta is `acos(2*k/(2*n+1))`, and gamma is `2*pi*k/phi`, where `phi=(1+sqrt(5))/2`. This gives `2*n+1` points.

For `'zcw'`, N is `fibonacci(n+2)`; k runs from 0 to N-1, beta is `acos(2*k/N-1)`, and gamma is `2*pi*k*fibonacci(n)/N`. For `'zcwn'`, N=n; the same beta rule uses n in place of N, and gamma is `2*pi*k/phi^2`. Both have N points. The formulas and ordering are fixed by the input; the function uses no random sampling.

The header help text lists `'fibonacci'`, but the executable switch label is `'fib'`; the literal `'fibonacci'` is not a switch case and reaches the unsupported-type error.

## Weights and plotting

The tessellation is calculated only when more than three outputs are requested or when the function is called with no output. Point coordinates are formed on the unit sphere and passed to `voronoisphere`; its body-angle weights are divided by `4*pi`, so the returned weights are normalised to the full sphere. A call with no outputs also sends the points and tessera to `grid_plot`. Requesting only the three angle outputs avoids this tessellation work.
