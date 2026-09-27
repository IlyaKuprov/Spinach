# examples/nmr_paramag/point_vs_distr.m

- Signature: `point_vs_distr()`

## Purpose

Comparison of point and multipole PCS fits when the electron probability density is spatially distributed: the synthetic density is a mixture of four randomly positioned Gaussians. The example samples nuclei in a spherical shell around the origin and fits the same PCS data with a point model and a multipole model through rank 2. The source cites the [model paper](http://dx.doi.org/10.1039/c6cp05437d). Calculation time: minutes.

## Physical / mathematical content

The four Gaussian centers are randomized within a 3 Å cube, with `sigma=0.5` Å. The simulated PCS data are evaluated by `kpcs` using the FFT method; inverse fits use `ippcs` for the point model and `ilpcs` for ranks `[0 1 2]`.

## Numerical / algorithmic content

The example uses a 64-point grid in each dimension and places 100 nuclei at random radii from 5 to 15 Å. It zero-pads the density by two volumes on each side before the FFT-based calculation, then reports fitted parameters, residuals, and multipole moments.

## Implementation structure

- Generate a four-Gaussian electron density and a random set of nuclear coordinates.
- Compute distributed PCS with the FFT solution and fit point and rank-2 multipole models.
- Compare fit residuals and recovered moments with the known input values.
