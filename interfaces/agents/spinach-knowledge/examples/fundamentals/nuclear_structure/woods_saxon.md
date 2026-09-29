# examples/fundamentals/nuclear_structure/woods_saxon.m

- MATLAB implementation: [examples/fundamentals/nuclear_structure/woods_saxon.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/nuclear_structure/woods_saxon.m)

- Signature: `woods_saxon(mass_number,level_number)`

## Purpose

A pedagogical finite-difference calculation of selected single-nucleon eigenstates in a three-dimensional central Woods-Saxon potential. The MATLAB source describes femtometres and MeV as its non-SI units and converts the kinetic prefactor from SI constants.

## Model and inputs

The radius is `r_nuc = r0 * mass_number^(1/3)`, with `r0 = 1.25 fm`. The radial potential is `V(R) = -V0 / (1 + exp((R-r_nuc)/a))`, with `V0 = 50 MeV` and `a = 0.5 fm`. If called with no arguments, the function sets `mass_number = 20` and `level_number = 5`; the source has no separate defaults for a partially supplied argument list and no explicit input validation.

The computational box spans `[-3*r_nuc, 3*r_nuc]` along each Cartesian axis, with 50 points per axis. The grid is constructed on the periodic grid used by `fdlap`; its last coordinate is one grid interval short of the upper box boundary. The kinetic operator is `L = premult * fdlap([50 50 50], box_sizes, 5)`, where `premult = 1e-6*1e30*hbar^2/(2*pmass*eV)`; the Hamiltonian is assembled as `H = diag(V(:)) - L`.

## Calculation and reported output

The code asks `eigs(H,level_number,'smallestreal')` for the requested eigenpairs, prints their energies in MeV, and makes volume plots of the potential and the real part of eigenvector `level_number`. Plot extents are adjusted to end at the last grid point. The source gives no reference energy values, convergence comparison, or pass/fail assertion, so the page does not treat the calculation as a validated spectrum.

## Scope and limitations

This is the single-particle central-potential model actually assembled in the source: a radial Woods-Saxon term plus a discretised kinetic term. It includes no explicit spin-orbit or other interaction term. Grid resolution, box size, and the finite-difference stencil are fixed in the example; no convergence or boundary-sensitivity study is reported. No MATLAB execution is claimed here.
