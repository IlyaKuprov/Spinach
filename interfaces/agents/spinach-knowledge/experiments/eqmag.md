# experiments/eqmag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/eqmag.m`
- Signature: `magn=eqmag(spin_system,parameters)`
- Total lines: 120

## Purpose

Computes the molar magnetization vector at the thermal equilibrium at the temperature specified in inter.temperature and magnetic field spe- cified in sys.magnet (assumed to be along the Z-axis), averaged over system orientations using the spherical grid specified. Syntax: magn=eqmag(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 32-33: Get the number of spins; implemented by `nspins=spin_system.comp.nspins`.
- Lines 35-36: Get the g-tensor for each spin in Bohr magneton units; implemented by `for n=nspins:-1:1`.
- Lines 40-41: Get Sx Sy Sz operators for each spin; implemented by `for n=nspins:-1:1`.
- Lines 47-48: Define the problem dimension; implemented by `problem_dim=prod(spin_system.comp.mults)`.
- Lines 50-51: Get Hamiltonian components; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 53-55: Load the spherical grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid '.mat'])`.
- Lines 59-60: Get the magnetization vector started; implemented by `magn=[0 0 0]`.
- Lines 62-63: Loop over the grid; implemented by `parfor k=1:numel(weights)`.
- Lines 65-66: Get the rotation matrix; implemented by `R=euler2dcm([alphas(k) betas(k) gammas(k)])`.
- Lines 68-69: Start magnetic moment operators; implemented by `mx=sparse(problem_dim,problem_dim)`.
- Lines 73-74: Compute the detection operator; implemented by `for n=1:nspins`.
- Lines 76-77: Rotate g-tensor (because we rotate the molecule, not the field); implemented by `g_rot=sparse(R*g{n}*R')`.
- Lines 79-80: Build magnetic moment operators; implemented by `mx=mx-g_rot(1,1)*Sx{n}-g_rot(1,2)*Sy{n}-g_rot(1,3)*Sz{n}`.
- Lines 86-87: Get the equilibrium density matrix; implemented by `rho=equilibrium(spin_system,I,Q,[alphas(k) betas(k) gammas(k)])`.
- Lines 89-90: Normalize density; implemented by `rho=rho./trace(rho)`.
- Lines 92-93: Get the magnetisation vector in Bohr magnetons per mol; implemented by `magn=magn+sph_grid.weights(k)*[trace(rho*mx) trace(rho*my) trace(rho*mz)]`.
- Lines 97-98: Keep the real part; implemented by `magn=real(magn)`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=nspins:-1:1`.
- Line 41: `for` loop over `n=nspins:-1:1`.
- Line 63: `parfor` loop over `k=1:numel(weights)`.
- Line 74: `for` loop over `n=1:nspins`.

### Key state/data transformations

- Lines 33: computes `nspins` using `nspins=spin_system.comp.nspins`.
- Lines 37: computes `g{n}` using `g{n}=gtensorof(spin_system,n)`.
- Lines 42: computes `Sx{n}` using `Sx{n}=operator(spin_system,{'Lx'},{n})`.
- Lines 43: computes `Sy{n}` using `Sy{n}=operator(spin_system,{'Ly'},{n})`.
- Lines 44: computes `Sz{n}` using `Sz{n}=operator(spin_system,{'Lz'},{n})`.
- Lines 48: computes `problem_dim` using `problem_dim=prod(spin_system.comp.mults)`.
- Lines 51: computes `[I,Q]` using `[I,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 54-55: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid '.mat'])`.
- Lines 56: computes `weights` using `weights=sph_grid.weights; alphas=sph_grid.alphas`.
- Lines 57: computes `betas` using `betas=sph_grid.betas; gammas=sph_grid.gammas`.
- Lines 60: computes `magn` using `magn=[0 0 0]`.
- Lines 66: computes `R` using `R=euler2dcm([alphas(k) betas(k) gammas(k)])`.
- Lines 69: computes `mx` using `mx=sparse(problem_dim,problem_dim)`.
- Lines 70: computes `my` using `my=sparse(problem_dim,problem_dim)`.
- Lines 71: computes `mz` using `mz=sparse(problem_dim,problem_dim)`.
- Lines 77: computes `g_rot` using `g_rot=sparse(R*g{n}*R')`.
- Lines 87: computes `rho` using `rho=equilibrium(spin_system,I,Q,[alphas(k) betas(k) gammas(k)])`.

### Local helper functions

- Line 103: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
  - Representative operation: `error('zeeman-hilb formalism is required.')`.

## Parameters / inputs

- parameters.grid -spherical grid for averaging

## Outputs

- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.

## Implementation structure

- Computes the molar magnetization vector at the thermal equilibrium at
- the temperature specified in inter.temperature and magnetic field spe-
- cified in sys.magnet (assumed to be along the Z-axis), averaged over
- system orientations using the spherical grid specified. Syntax:
- magn=eqmag(spin_system,parameters)
- parameters.grid -spherical grid for averaging
- magn -molar magnetization vector [Mx My Mz] in [Na*mu_bohr]
- Note: the use of bas.formalism='zeeman-hilb' is required.
- Note: Spinach uses NMR convention for the exchange coupling: exchange
- interaction term in the Hamiltonian is 2*pi*J*(LxSx+LySy+LzSz)
- where J is in Hz.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `load()`, `euler2dcm()`, `alphas()`, `betas()`, `gammas()`, `g_rot()`, `equilibrium()`, `strcmp()`, `isfield()`, `ischar()`.
