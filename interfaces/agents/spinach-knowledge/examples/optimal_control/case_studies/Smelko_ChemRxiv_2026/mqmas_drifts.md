# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m`
- Signature: `drifts=mqmas_drifts(spin_system,parameters)`
- Total lines: 187

## Purpose

Drift Hamiltonians of a quadrupolar nucleus under magic angle spin- ning, resolved in rotor phase, for every combination of a two-angle powder grid orientation and an initial rotor phase. The Hamiltonian is taken to second order in the rotating frame of the nucleus, and is held constant within each rotor phase tick. Syntax: drifts=mqmas_drifts(spin_system,parameters)

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 47-48: Get the laboratory frame Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 50-51: Get the carrier Hamiltonian; implemented by `C=carrier(spin_system,parameters.spins{1})`.
- Lines 53-55: Load the spherical integration grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas')`.
- Lines 57-60: Get rotor axis orientation; implemented by `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 63-64: Rotor phases at tick midpoints; implemented by `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 66-67: Tick shift between initial rotor phases; implemented by `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 69-70: Silence the workers; implemented by `spin_system.sys.output='hush'`.
- Lines 72-73: Parallel loop over grid orientations; implemented by `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 76-77: Preallocate the rotor stack; implemented by `stack=cell(1,parameters.n_ticks)`.
- Lines 79-80: Loop over rotor phase ticks; implemented by `for k=1:parameters.n_ticks`.
- Lines 82-83: Start with the isotropic part; implemented by `stack{k}=H`.
- Lines 85-86: Loop over spherical ranks; implemented by `for r=1:numel(Q)`.
- Lines 88-89: Compute crystallite orientation; implemented by `D_mol2rot=wigner(r,0,sph_grid.betas(n),sph_grid.alphas(n))`.
- Lines 91-92: Compute rotor axis tilt; implemented by `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 94-95: Compute rotor rotation; implemented by `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 97-98: Compose rotations; implemented by `D_comp=D_lab2rot*D_rotor*D_mol2rot`.
- Lines 100-101: Build the anisotropic part; implemented by `for p=1:(2*r+1)`.

### Control flow inferred from the code

- Line 74: `parfor` loop over `n=1:n_orients`.
- Line 80: `for` loop over `k=1:parameters.n_ticks`.
- Line 86: `for` loop over `r=1:numel(Q)`.
- Line 101: `for` loop over `p=1:(2*r+1)`.
- Line 102: `for` loop over `q=1:(2*r+1)`.
- Line 117: `for` loop over `j=1:parameters.n_phases`.

### Key state/data transformations

- Lines 48: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 51: computes `C` using `C=carrier(spin_system,parameters.spins{1})`.
- Lines 54-55: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas')`.
- Lines 58-60: computes `[rotor_phi,rotor_theta,~]` using `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 61: computes `rotor_theta` using `rotor_theta=pi/2-rotor_theta`.
- Lines 64: computes `rotor_phases` using `rotor_phases=2*pi*((1:parameters.n_ticks)-0.5)/parameters.n_ticks`.
- Lines 67: computes `tick_shift` using `tick_shift=parameters.n_ticks/parameters.n_phases`.
- Lines 70: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 73: computes `n_orients` using `n_orients=numel(sph_grid.alphas); orient_drifts=cell(1,n_orients)`.
- Lines 77: computes `stack` using `stack=cell(1,parameters.n_ticks)`.
- Lines 83: computes `stack{k}` using `stack{k}=H`.
- Lines 89: computes `D_mol2rot` using `D_mol2rot=wigner(r,0,sph_grid.betas(n),sph_grid.alphas(n))`.
- Lines 92: computes `D_lab2rot` using `D_lab2rot=wigner(r,rotor_phi,rotor_theta,0)`.
- Lines 95: computes `D_rotor` using `D_rotor=wigner(r,0,0,rotor_phases(k))`.
- Lines 98: computes `D_comp` using `D_comp=D_lab2rot*D_rotor*D_mol2rot`.
- Lines 116: computes `members` using `members=cell(1,parameters.n_phases)`.
- Lines 118-119: computes `tick_idx` using `tick_idx=mod((0:(parameters.n_slices-1))+(j-1)*tick_shift, parameters.n_ticks)+1`.
- Lines 120: computes `members{j}` using `members{j}=stack(tick_idx)`.

### Local helper functions

- Line 132: `grumble()` — `function grumble(spin_system,parameters)`.
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
- Note: the rotor turns the first Euler angle, therefore the azimuth
- of each grid point is placed into the third Euler angle and
- the first one is the initial rotor phase. The spin system must
- carry laboratory frame assumptions, call assume() first.

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

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `carrier()`, `load()`, `cart2sph()`, `wigner()`, `rotframe()`, `isfield()`, `iscell()`, `ischar()`, `ismember()`, `isrow()`, `isscalar()`.
