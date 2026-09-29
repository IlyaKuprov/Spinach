# experiments/pseudocon/lpcs.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/lpcs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=lpcs.m)

## Purpose

Computes theoretical pseudocontact shifts from multipole moments of the paramagnetic-centre probability density. It applies the paper’s Equation 33 under the stated far-field assumption that every nucleus lies outside the density’s bounding sphere; it is not a general near-field density solver.

## Inputs and moment convention

Call as `theo_pcs = lpcs(nxyz,mxyz,ranks,Ilm,chi)`. `nxyz` is a finite real N-by-3 array of nuclear coordinates in Å; `mxyz` is a finite real 1-by-3 paramagnetic-centre coordinate in Å. `ranks` is a vector of finite non-negative integer multipole ranks. `Ilm` is a cell array with one finite real vector of length `2*L+1` for each rank `L`.

For rank `L`, the supplied coefficients are spherical-harmonic moments: integrate the probability density times `Y_Lm` and `r^L` over volume. The source’s component ordering includes `Ilm = N/(2*sqrt(pi))` at rank zero; rank one is `[imag(I11), I10, real(I11)]`; rank two is `[imag(I22), imag(I21), I20, real(I21), real(I22)]`. Higher-rank vectors follow the corresponding real/imaginary component ordering.

`chi` may be a real 3-by-3 susceptibility tensor or its five independent elements, in Å³. The five-element form is expanded as a symmetric traceless tensor with the third diagonal equal to minus the sum of the first two. Coordinates are shifted by `mxyz` before conversion to spherical coordinates.

## Calculation and output

The routine combines each supplied multipole rank with the rank-2 anisotropic susceptibility components and spherical harmonics. A rank-`L` term has radial dependence `r^(-L-3)` and harmonic degree `L+2`; the summed real value is scaled by `1e6` to return `theo_pcs` in ppm at each nuclear coordinate.

## References

- [10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D)
- [Spin Dynamics Wiki: lpcs.m](https://spindynamics.org/wiki/index.php?title=lpcs.m)