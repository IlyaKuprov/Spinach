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
