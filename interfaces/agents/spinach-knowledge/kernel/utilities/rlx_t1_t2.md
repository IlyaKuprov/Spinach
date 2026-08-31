# kernel/utilities/rlx_t1_t2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_t1_t2.m`
- Signature: `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)`
- Total lines: 191

## Purpose

Extended T1/T2 relaxation model returning the relaxation super- operators separately for the longitudinal and the transverse states. Syntax: [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(spin_system)`.
- Lines 37-38: Compute ranks and projections; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 40-41: Preallocate relaxation rates; implemented by `r1_rates=zeros(size(spin_system.comp.isotopes))`.
- Lines 44-45: Fill in R1 relaxation rates; implemented by `for n=1:numel(spin_system.comp.isotopes)`.
- Lines 47-48: Get the specification; implemented by `current_r1_rate=spin_system.rlx.r1_rates{n}`.
- Lines 50-51: Scalar relaxation rate; implemented by `if isnumeric(current_r1_rate)&&isscalar(current_r1_rate)`.
- Lines 53-54: Simply assign the rate; implemented by `r1_rates(n)=current_r1_rate`.
- Lines 56-59: 3x3 tensor specification; implemented by `elseif isnumeric(current_r1_rate)&& (size(current_r1_rate,1)==3)&& (size(current_r1_rate,2)==3)`.
- Lines 61-62: Make sure the angles are specified; implemented by `if ~exist('euler_angles','var')`.
- Lines 66-67: Compute orientation ort (this matches alphas=0 of two-angle grids); implemented by `ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3))`.
- Lines 69-70: Get the rate at the current orientation; implemented by `r1_rates(n)=ort*current_r1_rate*ort'`.
- Lines 72-73: Function handle specification; implemented by `elseif isa(current_r1_rate,'function_handle')`.
- Lines 80-81: Call the function handle; implemented by `r1_rates(n)=current_r1_rate(euler_angles(1),euler_angles(2),euler_angles(3))`.
- Lines 85-86: Complain and bomb out; implemented by `error('unknown R1 rate specification.')`.
- Lines 92-93: Fill in R2 relaxation rates; implemented by `for n=1:numel(spin_system.comp.isotopes)`.
- Lines 95-96: Get the specification; implemented by `current_r2_rate=spin_system.rlx.r2_rates{n}`.
- Lines 98-99: Scalar relaxation rate; implemented by `if isnumeric(current_r2_rate)&&isscalar(current_r2_rate)`.
- Lines 101-102: Simply assign the rate; implemented by `r2_rates(n)=current_r2_rate`.

### Control flow inferred from the code

- Line 45: `for` loop over `n=1:numel(spin_system.comp.isotopes)`.
- Line 51: conditional branch on `isnumeric(current_r1_rate)&&isscalar(current_r1_rate)`.
- Line 62: conditional branch on `~exist('euler_angles','var')`.
- Line 76: conditional branch on `~exist('euler_angles','var')`.
- Line 93: `for` loop over `n=1:numel(spin_system.comp.isotopes)`.
- Line 99: conditional branch on `isnumeric(current_r2_rate)&&isscalar(current_r2_rate)`.
- Line 110: conditional branch on `~exist('euler_angles','var')`.
- Line 124: conditional branch on `~exist('euler_angles','var')`.
- Line 141: conditional branch on `any(~isreal(r1_rates),'all')||any(~isreal(r2_rates),'all')`.
- Line 151: `parfor` loop over `n=1:matrix_dim`.

### Key state/data transformations

- Lines 38: computes `[L,M]` using `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 41: computes `r1_rates` using `r1_rates=zeros(size(spin_system.comp.isotopes))`.
- Lines 42: computes `r2_rates` using `r2_rates=zeros(size(spin_system.comp.isotopes))`.
- Lines 48: computes `current_r1_rate` using `current_r1_rate=spin_system.rlx.r1_rates{n}`.
- Lines 54: computes `r1_rates(n)` using `r1_rates(n)=current_r1_rate`.
- Lines 67: computes `ort` using `ort=[0 0 1]*euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3))`.
- Lines 96: computes `current_r2_rate` using `current_r2_rate=spin_system.rlx.r2_rates{n}`.
- Lines 102: computes `r2_rates(n)` using `r2_rates(n)=current_r2_rate`.
- Lines 146: computes `matrix_dim` using `matrix_dim=size(spin_system.bas.basis,1)`.
- Lines 147: computes `r1_diagonal` using `r1_diagonal=zeros(matrix_dim,1)`.
- Lines 148: computes `r2_diagonal` using `r2_diagonal=zeros(matrix_dim,1)`.
- Lines 154: computes `local_r1_rates` using `local_r1_rates=r1_rates`.
- Lines 155: computes `local_r2_rates` using `local_r2_rates=r2_rates`.
- Lines 161: computes `r1_spins` using `r1_spins=(~logical(M(n,:)))&mask`.
- Lines 162: computes `r1_sum` using `r1_sum=sum(local_r1_rates(r1_spins))`.
- Lines 165: computes `r2_spins` using `r2_spins=(logical(M(n,:)))&mask`.
- Lines 166: computes `r2_sum` using `r2_sum=sum(local_r2_rates(r2_spins))`.
- Lines 169: computes `r1_diagonal(n)` using `r1_diagonal(n)=r1_sum`.

### Local helper functions

- Line 181: `grumble()` — `function grumble(spin_system)`. The state that separates its scholars from its warriors will have its thinking done by cowards and its fighting by fools.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
  - Representative operation: `error('this function is only available in sphten-liouv formalism.')`.

## Parameters / inputs

- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- red when R1 and/or R2 rates had been spe-
- cified as 3x3 tensor or a function handle,
- this argument has no effect for R1 and R2
- rates specified as scalars.

## Outputs

- R1Op -relaxation superoperator containing
- all longitudinal relaxation terms
- R2Op -relaxation superoperator containing
- all transverse relaxation terms
- Note: multi-spin orders relax at the sum of the rates of
- their constituent single-spin orders.

## Implementation structure

- Extended T1/T2 relaxation model returning the relaxation super-
- operators separately for the longitudinal and the transverse
- states. Syntax:
- [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)
- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- red when R1 and/or R2 rates had been spe-
- cified as 3x3 tensor or a function handle,
- this argument has no effect for R1 and R2
- rates specified as scalars.
- R1Op -relaxation superoperator containing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lin2lm()`, `isscalar()`, `r1_rates()`, `exist()`, `euler2dcm()`, `euler_angles()`, `current_r1_rate()`, `r2_rates()`, `current_r2_rate()`, `any()`, `logical()`, `local_r1_rates()`, `local_r2_rates()`, `r1_diagonal()`, `r2_diagonal()`.
