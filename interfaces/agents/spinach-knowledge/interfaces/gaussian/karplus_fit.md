# interfaces/gaussian/karplus_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gaussian/karplus_fit.m`
- Signature: `[A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)`
- Total lines: 124

## Purpose

Fits a Karplus curve to a Gaussian dihedral angle scan. Syntax: [A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- dir_path -path to the directory containing the
- Gaussian logs
- atoms -a cell array of 4-element vectors
- specifying atmos making up the dihe-
- dral angles of interest

## Outputs

- A,B,C -coefficients for A+B*cos(phi)+C*cos(phi)^2
- As,Bs,Vs -standard deviations of those coefficients
- The directory specified in the first argument should contain
- a series of Gaussian J-coupling calculation logs that differ
- only in the value of the dihedral angle in question.

## Implementation structure

- Fits a Karplus curve to a Gaussian dihedral angle scan. Syntax:
- [A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)
- dir_path -path to the directory containing the
- Gaussian logs
- atoms -a cell array of 4-element vectors
- specifying atmos making up the dihe-
- dral angles of interest
- A,B,C -coefficients for A+B*cos(phi)+C*cos(phi)^2
- As,Bs,Vs -standard deviations of those coefficients
- The directory specified in the first argument should contain
- a series of Gaussian J-coupling calculation logs that differ
- only in the value of the dihedral angle in question.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dir()`, `gparse()`, `logfiles()`, `phi()`, `dihedral()`, `isnan()`, `cosd()`, `result()`, `vec_res_sq()`, `jacobianest()`, `sum_res_sq()`, `inv()`, `stdevs()`, `kfigure()`, `kxlabel()`.
