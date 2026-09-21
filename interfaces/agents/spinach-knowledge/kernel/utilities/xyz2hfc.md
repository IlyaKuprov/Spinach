# kernel/utilities/xyz2hfc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/xyz2hfc.m`
- Signature: `A=xyz2hfc(exyz,nxyz,isotope)`
- Total lines: 80

## Purpose

Converts point electron and nuclear coordinates into a hyper- fine interaction tensor. Syntax: A=xyz2hfc(exyz,nxyz,isotope)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- exyz -Cartesian coordinates of the electron,
- a 1x3 row vector in Angstrom
- nxyz -Cartesian coordinates of the nucleus,
- a 1x3 row vector in Angstrom
- isotope -isitope specification, e.g. '13C'

## Outputs

- A -hyperfine coupling tensor, Gauss
- Note: Gauss units are used for hyperfine couplings because
- they do not depend on the electron g-tensor.
- Note: the tensor returned is the one that enters the spin
- Hamiltonian as S*A*I; it does not scale with the num-
- ber of unpaired electrons because the electron spin
- operator already carries that magnitude.

## Implementation structure

- Converts point electron and nuclear coordinates into a hyper-
- fine interaction tensor. Syntax:
- A=xyz2hfc(exyz,nxyz,isotope)
- exyz -Cartesian coordinates of the electron,
- a 1x3 row vector in Angstrom
- nxyz -Cartesian coordinates of the nucleus,
- isotope -isitope specification, e.g. '13C'
- A -hyperfine coupling tensor, Gauss
- Note: Gauss units are used for hyperfine couplings because
- they do not depend on the electron g-tensor.
- Note: the tensor returned is the one that enters the spin
- Hamiltonian as S*A*I; it does not scale with the num-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `isequal()`, `ischar()`.
