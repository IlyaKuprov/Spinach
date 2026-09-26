# kernel/plotting/efg_display.m

- Signature: `efg_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Electric field gradient tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors only): 1. A unit sphere in a Cartesian space is scaled by abs(Axx) in the x direction, abs(Ayy) in the y direction and abs(Azz) in the z direction, where Axx, Ayy, Azz are the eigenvalues of the CST ten- sor in units of ppm. 2. A set of axes is drawn inside the sphere with a red axis for a positive eig

## Physical / mathematical content

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Syntax

```matlab
efg_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- props -output of c2spinach or gparse
- atoms -a cell array of element symbols
- or a vector of integers, indica-
- ting the atoms for which EFG
- tensors should be visualised,
- e.g. {'N','O'} or [1 2 5]
- scaling -a factor to scale the tensors
- by for visualisation
- conmatrix -binary connectivity matrix, 1
- if a pair of atoms should be
- connected by a bond. If an em-
- pty vector is supplied, 1.6
- Angstrom cutoff distance is used
- options.style -'ellipsoids' or 'harmonics'
- options.kill_iso -set to true() to eliminate the
- isotropic parts of tensors be-
- fore plotting
- options.numbers -set to true() to display atom
- numbers
- options.symbols -set to false() to not display
- atom symbols

## Implementation structure

- Electric field gradient tensors and their eigensystems. Two
- styles are implemented:
- A. Ellipsoids (symmetric tensors only):
- 1. A unit sphere in a Cartesian space is scaled by
- abs(Axx) in the x direction, abs(Ayy) in the y
- direction and abs(Azz) in the z direction, where
- Axx, Ayy, Azz are the eigenvalues of the CST ten-
- sor in units of ppm.
- 2. A set of axes is drawn inside the sphere with a
- red axis for a positive eigenvalue, and a blue
- axis for a negative one.
- 3. The sphere is translated to the point of corres-
