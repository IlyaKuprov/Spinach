# examples/fundamentals/nuclear_structure/woods_saxon.m

- Signature: `woods_saxon(mass_number,level_number)`

## Purpose

Calculates selected single-nucleon eigenstates in a three-dimensional Woods–Saxon potential. Distances are in femtometres and energies in MeV, with SI constants used in the kinetic-energy conversion.

## Potential and discretisation

The constants are r0=1.25 fm, V0=50 MeV, and a=0.5 fm; the nuclear radius is r_nuc=r0*mass_number^(1/3). The potential on the grid is V=−V0/(1+exp((R-r_nuc)/a)). With no inputs, mass_number=20 and level_number=5. The periodic finite-difference Laplacian uses a 50-by-50-by-50 grid in a box spanning −3*r_nuc to +3*r_nuc on each axis and fdlap(box_npts,box_sizes,5). Its prefactor is 1e-6*1e30*hbar^2/(2*pmass*eV), with hbar=1.054571817e-34, pmass=1.6726219e-27, and eV=1.602176634e-19 in SI units.

## Eigenproblem and plots

The Hamiltonian is H=spdiags(V(:),0,prod(box_npts),prod(box_npts))-L. The code obtains level_number eigenpairs with eigs(H,level_number,'smallestreal') and displays their energies in MeV. It plots the potential and the real part of eigenvector level_number reshaped onto the grid; the plotted extents end at the last grid point, not the upper box boundary.
