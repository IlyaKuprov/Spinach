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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 61-62: Check consistency; implemented by `grumble(parameters)`.
- Lines 64-65: Get translation generators; implemented by `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`.
- Lines 67-68: Count voxels; implemented by `nvoxels=prod(parameters.npts)`.
- Lines 70-71: Default is no flow on X; implemented by `if ~isfield(parameters,'u'), parameters.u=0; end`.
- Lines 73-74: Flow in X direction; implemented by `if isscalar(parameters.u)`.
- Lines 76-77: Uniform flow on X; implemented by `F=parameters.u*Fx`.
- Lines 81-83: Non-uniform flow on X; implemented by `F=spdiags(Fx*parameters.u(:),0,nvoxels,nvoxels)+ spdiags(parameters.u(:),0,nvoxels,nvoxels)*Fx`.
- Lines 87-88: Uniform X diffusion; implemented by `if isfield(parameters,'diff')`.
- Lines 90-91: Diffusion type; implemented by `if isscalar(parameters.diff)`.
- Lines 93-94: Add isotropic diffusion on X; implemented by `F=F-1i*parameters.diff*Fx*Fx`.
- Lines 98-99: Add anisotropic diffusion on X; implemented by `F=F-1i*parameters.diff(1,1)*Fx*Fx`.
- Lines 105-106: Non-uniform X diffusion; implemented by `if isfield(parameters,'dxx')`.
- Lines 108-109: Add non-uniform diffusion on X; implemented by `F=F-1i*Fx*spdiags(parameters.dxx(:),0,nvoxels,nvoxels)*Fx`.
- Lines 113-114: When sample has Y dimension; implemented by `if numel(parameters.npts)>1`.
- Lines 116-117: Default is no flow on Y; implemented by `if ~isfield(parameters,'v'), parameters.v=0; end`.
- Lines 119-120: Flow in Y direction; implemented by `if isscalar(parameters.v)`.
- Lines 122-123: Add uniform flow on Y; implemented by `F=F+parameters.v*Fy`.
- Lines 127-129: Add non-uniform flow on Y; implemented by `F=F+spdiags(Fy*parameters.v(:),0,nvoxels,nvoxels)+ spdiags(parameters.v(:),0,nvoxels,nvoxels)*Fy`.

### Control flow inferred from the code

- Line 71: conditional branch on `~isfield(parameters,'u'), parameters.u=0; end`.
- Line 74: conditional branch on `isscalar(parameters.u)`.
- Line 88: conditional branch on `isfield(parameters,'diff')`.
- Line 91: conditional branch on `isscalar(parameters.diff)`.
- Line 106: conditional branch on `isfield(parameters,'dxx')`.
- Line 114: conditional branch on `numel(parameters.npts)>1`.
- Line 117: conditional branch on `~isfield(parameters,'v'), parameters.v=0; end`.
- Line 120: conditional branch on `isscalar(parameters.v)`.
- Line 134: conditional branch on `isfield(parameters,'diff')`.
- Line 137: conditional branch on `isscalar(parameters.diff)`.
- Line 154: conditional branch on `all(isfield(parameters,{'dxx','dxy','dyx','dyy'}))`.
- Line 166: conditional branch on `numel(parameters.npts)>2`.
- Line 169: conditional branch on `~isfield(parameters,'w'), parameters.w=0; end`.
- Line 172: conditional branch on `isscalar(parameters.w)`.

### Key state/data transformations

- Lines 65: computes `[Fx,Fy,Fz]` using `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`.
- Lines 68: computes `nvoxels` using `nvoxels=prod(parameters.npts)`.
- Lines 77: computes `F` using `F=parameters.u*Fx`.
- Lines 227: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1)`.

### Local helper functions

- Line 233: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'npts')`.
  - Representative operation: `error('the number of points in the spatial grid must be specified in parameters.npts field.')`.

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
