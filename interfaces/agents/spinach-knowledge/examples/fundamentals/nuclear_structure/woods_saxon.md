# examples/fundamentals/nuclear_structure/woods_saxon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/nuclear_structure/woods_saxon.m`
- Signature: `woods_saxon(mass_number,level_number)`
- Total lines: 58

## Purpose

A loose implementation of single-nucleon Hamiltonian eigenfunction calculation in the three-dimensional Woods-Saxon potential. Units, when not SI, are femtometres and MeV.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Fundamental constants; implemented by `r0=1.25; V0=50; a=0.5`.
- Lines 15-16: Defaults; implemented by `if nargin==0`.
- Lines 20-21: Nuclear radius; implemented by `r_nuc=r0*mass_number^(1/3)`.
- Lines 23-24: Simulation box dimensions; implemented by `box_extents=3*[-r_nuc,r_nuc,-r_nuc,r_nuc,-r_nuc,r_nuc]`.
- Lines 28-29: Laplacian part; implemented by `premult=1e-6*1e30*hbar^2/(2*pmass*eV)`.
- Lines 32-33: Potential part; implemented by `X=linspace(box_extents(1),box_extents(2),box_npts(1))`.
- Lines 39-40: Plot the potential; implemented by `kfigure(); volplot(V,box_extents)`.
- Lines 44-45: Assemble the Hamiltonian; implemented by `H=spdiags(V(:),0,prod(box_npts),prod(box_npts))-L`.
- Lines 47-48: Get the state; implemented by `[psi,E]=eigs(H,level_number,'smallestreal')`.
- Lines 51-52: Plot the state; implemented by `psi=reshape(real(psi(:,level_number)),box_npts)`.

### Control flow inferred from the code

- Line 16: conditional branch on `nargin==0`.

### Key state/data transformations

- Lines 10: computes `r0` using `r0=1.25; V0=50; a=0.5`.
- Lines 11: computes `hbar` using `hbar=1.054571817e-34`.
- Lines 12: computes `pmass` using `pmass=1.6726219e-27`.
- Lines 13: computes `eV` using `eV=1.602176634e-19`.
- Lines 17: computes `mass_number` using `mass_number=20; level_number=5`.
- Lines 21: computes `r_nuc` using `r_nuc=r0*mass_number^(1/3)`.
- Lines 24: computes `box_extents` using `box_extents=3*[-r_nuc,r_nuc,-r_nuc,r_nuc,-r_nuc,r_nuc]`.
- Lines 25: computes `box_sizes` using `box_sizes=6*[r_nuc r_nuc r_nuc]`.
- Lines 26: computes `box_npts` using `box_npts=[50 50 50]`.
- Lines 29: computes `premult` using `premult=1e-6*1e30*hbar^2/(2*pmass*eV)`.
- Lines 30: computes `L` using `L=premult*fdlap(box_npts,box_sizes,5)`.
- Lines 33: computes `X` using `X=linspace(box_extents(1),box_extents(2),box_npts(1))`.
- Lines 34: computes `Y` using `Y=linspace(box_extents(3),box_extents(4),box_npts(2))`.
- Lines 35: computes `Z` using `Z=linspace(box_extents(5),box_extents(6),box_npts(3))`.
- Lines 36: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(X,Y,Z); R=sqrt(X.^2+Y.^2+Z.^2)`.
- Lines 37: computes `V` using `V=-V0./(1+exp((R-r_nuc)/a))`.
- Lines 41: computes `ktitle(['Woods-Saxon potential, M` using `ktitle(['Woods-Saxon potential, M=' num2str(mass_number)])`.
- Lines 45: computes `H` using `H=spdiags(V(:),0,prod(box_npts),prod(box_npts))-L`.

## Implementation structure

- A loose implementation of single-nucleon Hamiltonian eigenfunction
- calculation in the three-dimensional Woods-Saxon potential. Units,
- when not SI, are femtometres and MeV.
- Fundamental constants
- Defaults
- Nuclear radius
- Simulation box dimensions
- Laplacian part
- Potential part
- Plot the potential
- Assemble the Hamiltonian
- Get the state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fdlap()`, `box_extents()`, `box_npts()`, `kfigure()`, `volplot()`, `ktitle()`, `num2str()`, `kxlabel()`, `kylabel()`, `kzlabel()`, `spdiags()`, `psi()`.
