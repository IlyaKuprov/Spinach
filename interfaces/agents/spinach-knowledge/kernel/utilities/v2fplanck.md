# kernel/utilities/v2fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/v2fplanck.m`
- Signature: `F=v2fplanck(spin_system,parameters)`
- Total lines: 428

## Purpose

Translates a stationary 3D velocity field and a diffusion tensor field into a Fokker-Planck evolution generator. Syntax: F=v2fplanck(spin_system,parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.u -X components of the velocity vectors
- for each voxel in the sample, m/s;
- a scalar specifies spatially uni-
- form flow along X
- parameters.v -Y components of the velocity vectors
- for each voxel in the sample, m/s;
- a scalar specifies spatially uni-
- form flow along Y
- parameters.w -Z components of the velocity vectors
- for each voxel in the sample, m/s;
- a scalar specifies spatially uni-
- form flow along Z
- parameters.diff -diffusion coefficient or 3x3 tensor, m^2/s
- for situations when this parameter is the
- same in every voxel
- parameters.dxx -Cartesian components of the diffusion
- parameters.dxy tensor for each voxel of the sample
- ...
- parameters.dzz
- parameters.dims -dimensions of the 3D box, meters
- parameters.npts -number of points in each dimension
- of the 3D box
- parameters.deriv -{'fourier'} uses Fourier diffe-
- rentiation matrices; {'period',n}
- requests n-point central finite-
- difference matrices with periodic
- boundary conditions

## Outputs

- F -spatial dynamics generator
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].
- Note: polyadic objects are returned, use inflate() to get the
- corresponding sparse matrix.

## Implementation structure

- Translates a stationary 3D velocity field and a diffusion tensor
- field into a Fokker-Planck evolution generator. Syntax:
- F=v2fplanck(spin_system,parameters)
- parameters.u -X components of the velocity vectors
- for each voxel in the sample, m/s;
- a scalar specifies spatially uni-
- form flow along X
- parameters.v -Y components of the velocity vectors
- form flow along Y
- parameters.w -Z components of the velocity vectors
- form flow along Z
- parameters.diff -diffusion coefficient or 3x3 tensor, m^2/s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hydrodynamics()`, `isfield()`, `isscalar()`, `spdiags()`, `all()`, `clean_up()`, `opium()`, `isrow()`, `any()`, `iscolumn()`, `isequal()`, `num2str()`, `getfield()`, `eps()`, `d_off()`.
