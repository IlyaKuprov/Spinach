# kernel/grids/grid_igloo.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_igloo.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=grid_igloo.m)

## Purpose

Construct the igloo spherical grid described in Appendix A.2 of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2014.05.009).

## Inputs and outputs

- `n_long` must be a positive real integer. Despite the parameter name and help description, the implementation uses it as the number of beta rings, including both poles.
- The outputs `alps`, `bets`, and `gams` are N-by-1 columns of Euler angles in radians; `alps` is zero. N is the total number of ring points. `whts` contains one normalised spherical-area weight per point, and `vorn` contains the Voronoi tessera vertex data.

## Grid rule

The beta values are evenly spaced from 0 through pi: `beta_k=k*pi/(n_long-1)` for k=0,...,n_long-1 when n_long is greater than one (the one-ring case is the single value 0). At each beta, the number of points is `m_k=max(1,floor(2*(n_long-1)*sin(beta_k)+0.5))`. Their gamma values are equally spaced around [0,2*pi), with spacing `2*pi/m_k`; the 2*pi endpoint is omitted. Thus the total count is the sum of the m_k values. This ring layout is deterministic and contains no random sampling.

## Weights and plotting

Voronoi weights are computed only when more than three outputs are requested or when there are no outputs. The Cartesian unit-sphere points are passed to `voronoisphere`, and its body-angle weights are divided by `4*pi` to normalise by the full sphere. A no-output call also plots the points and tessera using `grid_plot`.
