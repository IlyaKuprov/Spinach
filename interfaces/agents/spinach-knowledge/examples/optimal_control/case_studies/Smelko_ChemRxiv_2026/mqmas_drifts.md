# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m`
- Signature: `drifts=mqmas_drifts(spin_system,parameters)`
- Total lines: 193

## Purpose

Drift Hamiltonians of a quadrupolar nucleus under magic angle spinning, resolved in rotor phase, for every combination of a two-angle powder grid orientation and an initial rotor phase. The Hamiltonian is taken to second order in the rotating frame of the nucleus, and is held constant within each rotor phase tick. Syntax: drifts=mqmas_drifts(spin_system,parameters) Parameters: parameters.spins - the nucleus, a cell array with one isotope string, e.g. {'27Al'} parameters.axis - spinning axis, a normalised row vector with three elements parameters.grid - two-angle powder grid name; the grid must have uniform weights because the ensemble average in optimcon is unweighted parameters.n_ticks - rotor phase ticks per rotor period parameters.n_phases - number of initial rotor phases, must be a divisor of parameters.n_ticks parameters.n_slices - number of ticks in the pulse Outputs: drifts - cell array over the ensemble, grid orientations in the outer index and initial rotor phases in the inner index; each element is a cell array of n_slices Hamiltonian matrices, one per tick of the pulse Note: the crystallite orientation uses all three Euler angles of the grid, as in singlerot.m; two-angle grids keep the azimuth in the third angle and have zero first angles, which the rotor phase then supplies. The spin system must carry laboratory frame assumptions, call assume() first.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 48-49: Get the laboratory frame Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 51-52: Get the carrier Hamiltonian; implemented by `C=carrier(spin_system,parameters.spins{1})`.
- Lines 54-56: Load the spherical integration grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 58-59: Uniformly weighted grids only, the ensemble average is unweighted; implemented by `if any(abs(sph_grid.weights-sph_grid.weights(1))>1e-12)`.
- Lines 63-66: Get rotor axis orientation; implemented by `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 69-70: Rotor phases at tick midpoints; implemented by `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 72-73: Tick shift between initial rotor phases; implemented by `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 75-76: Silence the workers; implemented by `spin_system.sys.output='hush'`.
- Lines 78-79: Parallel loop over grid orientations; implemented by `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 82-83: Preallocate the rotor stack; implemented by `stack=cell(1,parameters.n_ticks)`.
- Lines 85-86: Loop over rotor phase ticks; implemented by `for k=1:parameters.n_ticks`.
- Lines 88-89: Start with the isotropic part; implemented by `stack{k}=H`.
- Lines 91-92: Loop over spherical ranks; implemented by `for r=1:numel(Q)`.
- Lines 94-95: Compute crystallite orientation; implemented by `D_mol2rot=wigner(r,sph_grid.alphas(n),sph_grid.betas(n),sph_grid.gammas(n))`.
- Lines 97-98: Compute rotor axis tilt; implemented by `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 100-101: Compute rotor rotation; implemented by `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 103-104: Compose rotations; implemented by `D_comp=D_lab2rot*D_rotor*D_mol2rot`.

### Control flow inferred from the code

- Line 59: conditional branch on `any(abs(sph_grid.weights-sph_grid.weights(1))>1e-12)`.
- Line 80: `parfor` loop over `n=1:n_orients`.
- Line 86: `for` loop over `k=1:parameters.n_ticks`.
- Line 92: `for` loop over `r=1:numel(Q)`.
- Line 107: `for` loop over `p=1:(2*r+1)`.
- Line 108: `for` loop over `q=1:(2*r+1)`.
- Line 123: `for` loop over `j=1:parameters.n_phases`.

### Key state/data transformations

- Lines 49: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 52: computes `C` using `C=carrier(spin_system,parameters.spins{1})`.
- Lines 55-56: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 64-66: computes `[rotor_phi,rotor_theta,~]` using `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 67: computes `rotor_theta` using `rotor_theta=pi/2-rotor_theta`.
- Lines 70: computes `rotor_phases` using `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 73: computes `tick_shift` using `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 76: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 79: computes `n_orients` using `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 83: computes `stack` using `stack=cell(1,parameters.n_ticks)`.
- Lines 89: computes `stack{k}` using `stack{k}=H`.
- Lines 95: computes `D_mol2rot` using `D_mol2rot=wigner(r,sph_grid.alphas(n),sph_grid.betas(n),sph_grid.gammas(n))`.
- Lines 98: computes `D_lab2rot` using `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 101: computes `D_rotor` using `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 104: computes `D_comp` using `D_comp=D_lab2rot*D_rotor*D_mol2rot`.
- Lines 122: computes `members` using `members=cell(1,parameters.n_phases)`.
- Lines 124-125: computes `tick_idx` using `tick_idx=mod((0:(parameters.n_slices-1))+(j-1)*tick_shift, parameters.n_ticks)+1`.
- Lines 126: computes `members{j}` using `members{j}=stack(tick_idx)`.

### Local helper functions

- Line 138: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'spins')`.
  - Representative operation: `error('the nucleus must be specified in parameters.spins field.')`.

## Parameters / inputs

- parameters.spins -the nucleus, a cell array with one
- isotope string, e.g. {'27Al'}
- parameters.axis -spinning axis, a normalised row
- vector with three elements
- parameters.grid -two-angle powder grid name; the grid
- must have uniform weights because the
- ensemble average in optimcon is unweighted
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
- parameters.grid -two-angle powder grid name; the grid
- must have uniform weights because the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `carrier()`, `load()`, `any()`, `cart2sph()`, `wigner()`, `rotframe()`, `isfield()`, `iscell()`, `ischar()`, `ismember()`, `isrow()`.
