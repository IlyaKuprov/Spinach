# kernel/plotting/cst_display.m

- Signature: `cst_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Draws shielding tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors only): 1. A unit sphere in a Cartesian space is scaled by abs(Axx) in the x direction, abs(Ayy) in the y direction and abs(Azz) in the z direction, where Axx, Ayy, Azz are the eigenvalues of the CST ten- sor in units of ppm. 2. A set of axes is drawn inside the sphere with a red axis for a positive eigenvalue,

## Physical / mathematical content

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Syntax

```matlab
cst_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- props -output of gparse function
- atoms -a cell array of element symbols
- or a vector of integers, indica-
- ting the atoms for which shiel-
- ding tensors should be visuali-
- sed, e.g. {'C','H'} or [1 2 5]
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

- Draws shielding tensors and their eigensystems. Two styles are
- implemented:
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
