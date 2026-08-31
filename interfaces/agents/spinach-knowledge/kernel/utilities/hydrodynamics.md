# kernel/utilities/hydrodynamics.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/hydrodynamics.m`
- Signature: `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`
- Total lines: 165

## Purpose

A basic hydrodynamics infrastructure provider, returns first derivative operators with respect to the three sample coordi- nates. Syntax: [Fx,Fy,Fz]=hydrodynamics(parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(parameters)`.
- Lines 44-45: Build derivative operators; implemented by `switch parameters.deriv{1}`.
- Lines 49-50: Finite-difference derivatives; implemented by `if isscalar(parameters.npts)`.
- Lines 65-66: Fourier derivatives; implemented by `if isscalar(parameters.npts)`.
- Lines 81-82: Complain and bomb out; implemented by `error('unrecognized derivative operator type.')`.
- Lines 86-87: Kron up derivative operators; implemented by `Fx=[]; Fy=[]; Fz=[]`.
- Lines 101-102: Inflate derivative operators; implemented by `if ~ismember('polyadic',spin_system.sys.enable)`.

### Control flow inferred from the code

- Line 45: dispatches on `parameters.deriv{1}`; cases `'period'`, `'fourier'`.
- Line 50: conditional branch on `isscalar(parameters.npts)`.
- Line 53: conditional branch on `numel(parameters.npts)==2`.
- Line 57: conditional branch on `numel(parameters.npts)==3`.
- Line 66: conditional branch on `isscalar(parameters.npts)`.
- Line 69: conditional branch on `numel(parameters.npts)==2`.
- Line 73: conditional branch on `numel(parameters.npts)==3`.
- Line 88: conditional branch on `isscalar(parameters.npts)`.
- Line 91: conditional branch on `numel(parameters.npts)==2`.
- Line 95: conditional branch on `numel(parameters.npts)==3`.
- Line 102: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.

### Key state/data transformations

- Lines 51: computes `Dx` using `Dx=fdmat(parameters.npts(1),parameters.deriv{2},1)/(parameters.dims(1)/parameters.npts(1))`.
- Lines 55: computes `Dy` using `Dy=fdmat(parameters.npts(2),parameters.deriv{2},1)/(parameters.dims(2)/parameters.npts(2))`.
- Lines 60: computes `Dz` using `Dz=fdmat(parameters.npts(3),parameters.deriv{2},1)/(parameters.dims(3)/parameters.npts(3))`.
- Lines 67: computes `[~,Dx]` using `[~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx`.
- Lines 71: computes `[~,Dy]` using `[~,Dy]=fourdif(parameters.npts(2),1); Dy=(2*pi/parameters.dims(2))*Dy`.
- Lines 76: computes `[~,Dz]` using `[~,Dz]=fourdif(parameters.npts(3),1); Dz=(2*pi/parameters.dims(3))*Dz`.
- Lines 87: computes `Fx` using `Fx=[]; Fy=[]; Fz=[]`.
- Lines 93: computes `Fy` using `Fy=-1i*polyadic({{Dy,opium(parameters.npts(1),1)}})`.
- Lines 98: computes `Fz` using `Fz=-1i*polyadic({{Dz,opium(parameters.npts(2),1),opium(parameters.npts(1),1)}})`.

### Local helper functions

- Line 109: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'dims')`.
  - Representative operation: `error('sample dimensions must be specified in parameters.dims variable.')`.

## Parameters / inputs

- parameters.dims -dimensions of the sample (meters),
- one, two, or three-element row
- vector
- parameters.npts -number of points in each dimension
- of the sample, one, two, or three-
- element row vector
- parameters.deriv -{'fourier'} requests Fourier diffe-
- rentiation matrices; {'period',n}
- requests n-point central finite-
- difference matrices with periodic
- boundary conditions

## Outputs

- Fx, Fy, Fz -derivative matrices, SI units
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].
- Note: polyadic objects are returned, use inflate() to get the
- corresponding sparse matrix.

## Implementation structure

- A basic hydrodynamics infrastructure provider, returns first
- derivative operators with respect to the three sample coordi-
- nates. Syntax:
- [Fx,Fy,Fz]=hydrodynamics(parameters)
- parameters.dims -dimensions of the sample (meters),
- one, two, or three-element row
- vector
- parameters.npts -number of points in each dimension
- of the sample, one, two, or three-
- element row vector
- parameters.deriv -{'fourier'} requests Fourier diffe-
- rentiation matrices; {'period',n}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `fdmat()`, `fourdif()`, `polyadic()`, `opium()`, `ismember()`, `inflate()`, `isfield()`, `any()`, `iscell()`, `ischar()`, `strcmp()`.
