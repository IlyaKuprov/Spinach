# kernel/utilities/apodisation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/apodisation.m`
- Signature: `fid=apodisation(spin_system,fid,winfuns,fp_half)`
- Total lines: 253

## Purpose

Performs free induction decay apodisation. Supports free induction decays of any dimension. To satisfy Fourier transform symmetry requirements, the first elements of the FID in each dimension are divided by 2, except for singleton dimensions and those the user designates inactive. Syntax: fid=apodisation(spin_system,fid,winfuns)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `true()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- fid -the free induction decay. The function expects a column
- vector in the case of 1D FID, a 2D matrix with the time
- origin located at the (1,1) corner point in the case of
- a 2D FID, a 3D matrix with the time origin located at
- the (1,1,1) corner point in the case of a 3D FID, etc.
- winfuns -a cell array of window function specifications for each
- dimension of the FID in the format {{spec},{spec},...},
- omitting singleton dimensions. The following specifica-
- tions are supported:
- {} -do nothing in this dimension
- {'none'} -no window function, but divide the first
- point by 2 to satisfy the Fourier trans-
- form symmetry requirement
- {'crisp'} -multiplied by cos(x)^8 half-bell. First
- point has x=0, last point has x=pi/2.
- {'exp',k} -multiplied by exp(-k*x). First point has
- x=0, last point has x=1.
- {'gauss',k} -multiplied by exp(-k*(x.^2)). First point
- has x=0, last point has x=1.
- {'cos'} -multiplied by cos(x) half-bell. First po-
- int has x=0, last point has x=pi/2.
- {'sin'} -multiplied by sin(x) full bell. First po-
- int has x=0, last point has x=pi.
- {'sqcos'} -multiplied by cos(x).^2 half-bell. First
- point has x=0, last point has x=pi/2.
- {'sqsin'} -multiplied by sin(x).^2 full bell. First
- point has x=0, last point has x=pi.
- {'kaiser',k} -multiplied by a Kaiser function with the
- side lobe attenuation factor k. The peak
- of the Kaiser function is in the middle
- of the FID.
- {'bad-z1',k} -emulation of a misset Z1 shim, k is a di-
- mensionless constant proportional to the
- shim current, a good guess for 1H NMR at
- 600 MHz is 10.
- {'bad-z2',k} -emulation of a misset Z2 shim, k is a di-
- mensionless constant proportional to the
- shim current, a good guess for 1H NMR at
- 600 MHz is 40.
- fp_half -set to false() to disable dividing of the first points
- by 2, this is needed when multiple window functions are
- applied to the same dimension one after another

## Outputs

- fid -apodised free induction decay

## Implementation structure

- Performs free induction decay apodisation. Supports free induction decays
- of any dimension. To satisfy Fourier transform symmetry requirements, the
- first elements of the FID in each dimension are divided by 2, except for
- singleton dimensions and those the user designates inactive. Syntax:
- fid=apodisation(spin_system,fid,winfuns)
- fid -the free induction decay. The function expects a column
- vector in the case of 1D FID, a 2D matrix with the time
- origin located at the (1,1) corner point in the case of
- a 2D FID, a 3D matrix with the time origin located at
- the (1,1,1) corner point in the case of a 3D FID, etc.
- winfuns -a cell array of window function specifications for each
- dimension of the FID in the format {{spec},{spec},...},

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `exist()`, `true()`, `ndims()`, `rel_dims()`, `false()`, `cellfun()`, `setdiff()`, `fid()`, `report()`, `num2str()`, `kaiser()`, `transpose()`, `sinc()`, `fresnelc()`, `fresnels()`.
