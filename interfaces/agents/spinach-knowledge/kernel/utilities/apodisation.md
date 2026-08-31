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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 77-78: Check consistency; implemented by `grumble(fid,winfuns)`.
- Lines 80-81: Default is to halve first points; implemented by `if ~exist('fp_half','var'), fp_half=true(); end`.
- Lines 83-84: Find non-singleton dimensions; implemented by `rel_dims=true(1,ndims(fid))`.
- Lines 88-89: Exclude inactive dimensions; implemented by `inact_dims=find(cellfun(@isempty,winfuns))`.
- Lines 92-93: Factors of 2; implemented by `if fp_half`.
- Lines 96-97: Index first points along each dimension; implemented by `idx=repmat({':'},1,ndims(fid)); idx{dim}=1`.
- Lines 99-100: FFT symmetry requirement; implemented by `fid(idx{:})=fid(idx{:})/2`.
- Lines 102-104: Report to the user; implemented by `report(spin_system,['FID dimension ' num2str(dim) ', all first points divided by 2'])`.
- Lines 109-110: Window functions; implemented by `for n=1:numel(rel_dims)`.
- Lines 112-113: Get dimension and its point count; implemented by `dim=rel_dims(n); npts=size(fid,dim)`.
- Lines 115-116: Build window function; implemented by `switch winfuns{dim}{1}`.
- Lines 166-167: For all ye lazy PhD students out there; implemented by `x=linspace(0,1,npts); x=transpose(x(2:end))`.
- Lines 170-171: Avoid the singularity; implemented by `wf=[1; wf]`.
- Lines 175-177: Well burn my papers and call me a teaching fellow: turns out we need this too, you sloppy muppets; implemented by `x=linspace(0,1,npts); x=transpose(x(2:end))`.
- Lines 180-181: Symbolic toolbox is slow; implemented by `parfor k=1:numel(x)`.
- Lines 190-192: Complain and bomb out; implemented by `error(['window function type ' winfuns{dim}{1} ' is not implemented.'])`.
- Lines 196-197: Apply in the active dimension; implemented by `idx=ones(1,ndims(fid)); idx(dim)=npts`.
- Lines 200-202: Report to the user; implemented by `report(spin_system,['FID dimension ' num2str(dim) ', ' winfuns{dim}{1} ' window function applied.'])`.

### Control flow inferred from the code

- Line 81: conditional branch on `~exist('fp_half','var'), fp_half=true(); end`.
- Line 93: conditional branch on `fp_half`.
- Line 94: `for` loop over `dim=rel_dims`.
- Line 110: `for` loop over `n=1:numel(rel_dims)`.
- Line 116: dispatches on `winfuns{dim}{1}`; cases `'none'`, `'crisp'`, `'exp'`, `'gauss'`, `'cos'`, `'sin'`, `'sqcos'`, `'sqsin'`, `'kaiser'`, `'bad-z1'`, ….
- Line 181: `parfor` loop over `k=1:numel(x)`.

### Key state/data transformations

- Lines 84: computes `rel_dims` using `rel_dims=true(1,ndims(fid))`.
- Lines 85: computes `rel_dims(size(fid)<2)` using `rel_dims(size(fid)<2)=false()`.
- Lines 89: computes `inact_dims` using `inact_dims=find(cellfun(@isempty,winfuns))`.
- Lines 97: computes `idx` using `idx=repmat({':'},1,ndims(fid)); idx{dim}=1`.
- Lines 100: computes `fid(idx{:})` using `fid(idx{:})=fid(idx{:})/2`.
- Lines 113: computes `dim` using `dim=rel_dims(n); npts=size(fid,dim)`.
- Lines 120: computes `wf` using `wf=ones(npts,1)`.
- Lines 124: computes `x` using `x=linspace(0,pi/2,npts)`.
- Lines 130: computes `k` using `k=winfuns{dim}{2}`.
- Lines 182: computes `wf(k)` using `wf(k)=(fresnelc(x(k))+1i*fresnels(x(k)))/x(k)`.
- Lines 198: computes `fid` using `fid=fid.*reshape(wf,idx)`.

### Local helper functions

- Line 209: `grumble()` — `function grumble(fid,winfuns)`.
  - Representative operation: `if ~isnumeric(fid)`.
  - Representative operation: `error('fid must be a numeric array.')`.

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
