# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m`
- Signature: `drifts=mqmas_drifts(spin_system,parameters)`
- Total lines: 186

## Purpose

Drift Hamiltonians of a quadrupolar nucleus under magic angle spinning, resolved in rotor phase, for every combination of a two-angle powder grid orientation and an initial rotor phase. The Hamiltonian is taken to second order in the rotating frame of the nucleus, and is held constant within each rotor phase tick. Syntax: drifts=mqmas_drifts(spin_system,parameters) Parameters: parameters.spins - the nucleus, a cell array with one isotope string, e.g. {'27Al'} parameters.axis - spinning axis, a normalised row vector with three elements parameters.grid - two-angle powder grid name parameters.n_ticks - rotor phase ticks per rotor period parameters.n_phases - number of initial rotor phases, must be a divisor of parameters.n_ticks parameters.n_slices - number of ticks in the pulse Outputs: drifts - cell array over the ensemble, grid orientations in the outer index and initial rotor phases in the inner index; each element is a cell array of n_slices Hamiltonian matrices, one per tick of the pulse Note: the crystallite orientation uses all three Euler angles of the grid, as in singlerot.m; two-angle grids keep the azimuth in the third angle and have zero first angles, which the rotor phase then supplies. The spin system must carry laboratory frame assumptions, call assume() first.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 46-47: Get the laboratory frame Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 49-50: Get the carrier Hamiltonian; implemented by `C=carrier(spin_system,parameters.spins{1})`.
- Lines 52-54: Load the spherical integration grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas')`.
- Lines 56-59: Get rotor axis orientation; implemented by `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 62-63: Rotor phases at tick midpoints; implemented by `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 65-66: Tick shift between initial rotor phases; implemented by `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 68-69: Silence the workers; implemented by `spin_system.sys.output='hush'`.
- Lines 71-72: Parallel loop over grid orientations; implemented by `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 75-76: Preallocate the rotor stack; implemented by `stack=cell(1,parameters.n_ticks)`.
- Lines 78-79: Loop over rotor phase ticks; implemented by `for k=1:parameters.n_ticks`.
- Lines 81-82: Start with the isotropic part; implemented by `stack{k}=H`.
- Lines 84-85: Loop over spherical ranks; implemented by `for r=1:numel(Q)`.
- Lines 87-88: Compute crystallite orientation; implemented by `D_mol2rot=wigner(r,sph_grid.alphas(n),sph_grid.betas(n),sph_grid.gammas(n))`.
- Lines 90-91: Compute rotor axis tilt; implemented by `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 93-94: Compute rotor rotation; implemented by `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 96-97: Compose rotations; implemented by `D_comp=D_lab2rot*D_rotor*D_mol2rot`.
- Lines 99-100: Build the anisotropic part; implemented by `for p=1:(2*r+1)`.

### Control flow inferred from the code

- Line 73: `parfor` loop over `n=1:n_orients`.
- Line 79: `for` loop over `k=1:parameters.n_ticks`.
- Line 85: `for` loop over `r=1:numel(Q)`.
- Line 100: `for` loop over `p=1:(2*r+1)`.
- Line 101: `for` loop over `q=1:(2*r+1)`.
- Line 116: `for` loop over `j=1:parameters.n_phases`.

### Key state/data transformations

- Lines 47: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 50: computes `C` using `C=carrier(spin_system,parameters.spins{1})`.
- Lines 53-54: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas')`.
- Lines 57-59: computes `[rotor_phi,rotor_theta,~]` using `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 60: computes `rotor_theta` using `rotor_theta=pi/2-rotor_theta`.
- Lines 63: computes `rotor_phases` using `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 66: computes `tick_shift` using `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 69: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 72: computes `n_orients` using `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 76: computes `stack` using `stack=cell(1,parameters.n_ticks)`.
- Lines 82: computes `stack{k}` using `stack{k}=H`.
- Lines 88: computes `D_mol2rot` using `D_mol2rot=wigner(r,sph_grid.alphas(n),sph_grid.betas(n),sph_grid.gammas(n))`.
- Lines 91: computes `D_lab2rot` using `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 94: computes `D_rotor` using `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 97: computes `D_comp` using `D_comp=D_lab2rot*D_rotor*D_mol2rot`.
- Lines 115: computes `members` using `members=cell(1,parameters.n_phases)`.
- Lines 117-118: computes `tick_idx` using `tick_idx=mod((0:(parameters.n_slices-1))+(j-1)*tick_shift, parameters.n_ticks)+1`.
- Lines 119: computes `members{j}` using `members{j}=stack(tick_idx)`.

### Local helper functions

- Line 131: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'spins')`.
  - Representative operation: `error('the nucleus must be specified in parameters.spins field.')`.

## Parameters / inputs

- parameters.spins -the nucleus, a cell array with one
- isotope string, e.g. {'27Al'}
- parameters.axis -spinning axis, a normalised row
- vector with three elements
- parameters.grid -two-angle powder grid name
- parameters.n_ticks -rotor phase ticks per rotor period
- parameters.n_phases -number of initial rotor phases, must
- be a divisor of parameters.n_ticks
- parameters.n_slices -number of ticks in the pulse

## Outputs

- drifts -cell array over the ensemble, grid orientations in
- the outer index and initial rotor phases in the in-
- ner index; each element is a cell array of n_slices
- Hamiltonian matrices, one per tick of the pulse
- Note: the crystallite orientation uses all three Euler angles of the
- grid, as in singlerot.m; two-angle grids keep the azimuth in
- the third angle and have zero first angles, which the rotor
- phase then supplies. The spin system must carry laboratory
- frame assumptions, call assume() first.

## Implementation structure

- Drift Hamiltonians of a quadrupolar nucleus under magic angle spin-
- ning, resolved in rotor phase, for every combination of a two-angle
- powder grid orientation and an initial rotor phase. The Hamiltonian
- is taken to second order in the rotating frame of the nucleus, and
- is held constant within each rotor phase tick. Syntax:
- drifts=mqmas_drifts(spin_system,parameters)
- parameters.spins -the nucleus, a cell array with one
- isotope string, e.g. {'27Al'}
- parameters.axis -spinning axis, a normalised row
- vector with three elements
- parameters.grid -two-angle powder grid name
- parameters.n_ticks -rotor phase ticks per rotor period

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `carrier()`, `load()`, `cart2sph()`, `wigner()`, `rotframe()`, `isfield()`, `iscell()`, `ischar()`, `ismember()`, `isrow()`, `any()`.
