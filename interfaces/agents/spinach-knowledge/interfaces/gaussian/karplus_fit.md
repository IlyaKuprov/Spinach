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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(dir_path,atoms)`.
- Lines 33-34: Get all log files in the directory; implemented by `logfiles=dir([dir_path filesep '*.log'])`.
- Lines 36-37: Preallocate the arrays; implemented by `phi=zeros(numel(atoms),numel(logfiles))`.
- Lines 40-41: Parse the files; implemented by `parfor n=1:numel(logfiles)`.
- Lines 45-46: Extract parameters; implemented by `for n=1:numel(logfiles)`.
- Lines 61-62: Eliminate NaNs; implemented by `mask=isnan(phi)|isnan(J)`.
- Lines 65-66: Stretch the arrays; implemented by `phi=phi(:)'; J=J(:)'`.
- Lines 68-69: Rotate phi into [0,360]; implemented by `phi=mod(phi,360)`.
- Lines 71-72: Compute basis cosines; implemented by `cosine_2=cosd(phi).^2`.
- Lines 76-77: Run linear least squares; implemented by `result=[cosine_2' cosine_1' cosine_0']\J'`.
- Lines 80-81: Vector of squared residuals; implemented by `vec_res_sq=@(x)(J'-x(1)*cosine_2'-x(2)*cosine_1'-x(3)*cosine_0').^2`.
- Lines 83-84: Sum of squared residuals; implemented by `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 86-87: Compute the Jacobian at the optimal point; implemented by `jac=jacobianest(vec_res_sq,result)`.
- Lines 89-90: Get the Studentized residual; implemented by `sdr=sqrt(sum_res_sq(result)/(numel(J)-3))`.
- Lines 92-93: Get the standard deviations; implemented by `stdevs=sqrt(diag((sdr^2)*inv(jac'*jac))); %#ok<MINV>`.
- Lines 96-97: Plot Karplus curve; implemented by `kfigure(); plot(phi,J,'ro'); kgrid`.

### Control flow inferred from the code

- Line 41: `parfor` loop over `n=1:numel(logfiles)`.
- Line 46: `for` loop over `n=1:numel(logfiles)`.
- Line 48: `for` loop over `k=1:numel(atoms)`.
- Line 55: `for` loop over `k=1:numel(atoms)`.

### Key state/data transformations

- Lines 34: computes `logfiles` using `logfiles=dir([dir_path filesep '*.log'])`.
- Lines 37: computes `phi` using `phi=zeros(numel(atoms),numel(logfiles))`.
- Lines 38: computes `J` using `J=zeros(numel(atoms),numel(logfiles))`.
- Lines 42: computes `props{n}` using `props{n}=gparse([dir_path filesep logfiles(n).name])`.
- Lines 49-50: computes `phi(k,n)` using `phi(k,n)=dihedral(props{n}.std_geom(atoms{k}(1),:),props{n}.std_geom(atoms{k}(2),:), props{n}.std_geom(atoms{k}(3),:),props{n}.std_geom(atoms{k}(4),:))`.
- Lines 51: computes `J(k,n)` using `J(k,n)=props{n}.j_couplings(atoms{k}(1),atoms{k}(4))`.
- Lines 62: computes `mask` using `mask=isnan(phi)|isnan(J)`.
- Lines 63: computes `phi(mask)` using `phi(mask)=[]; J(mask)=[]`.
- Lines 72: computes `cosine_2` using `cosine_2=cosd(phi).^2`.
- Lines 73: computes `cosine_1` using `cosine_1=cosd(phi)`.
- Lines 74: computes `cosine_0` using `cosine_0=ones(size(phi))`.
- Lines 77: computes `result` using `result=[cosine_2' cosine_1' cosine_0']\J'`.
- Lines 78: computes `A` using `A=result(1); B=result(2); C=result(3)`.
- Lines 81: computes `vec_res_sq` using `vec_res_sq=@(x)(J'-x(1)*cosine_2'-x(2)*cosine_1'-x(3)*cosine_0').^2`.
- Lines 84: computes `sum_res_sq` using `sum_res_sq=@(x)sum(vec_res_sq(x))`.
- Lines 87: computes `jac` using `jac=jacobianest(vec_res_sq,result)`.
- Lines 90: computes `sdr` using `sdr=sqrt(sum_res_sq(result)/(numel(J)-3))`.
- Lines 93: computes `stdevs` using `stdevs=sqrt(diag((sdr^2)*inv(jac'*jac))); %#ok<MINV>`.

### Local helper functions

- Line 106: `grumble()` — `function grumble(dir_path,atoms)`. If men learn this, it will implant forgetfulness in their souls; they will cease to exercise memory because they rely on that which is written [...]
  - Representative operation: `if ~ischar(dir_path)`.
  - Representative operation: `error('dir_path must be a character string.')`.

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
