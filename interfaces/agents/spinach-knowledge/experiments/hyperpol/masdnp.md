# experiments/hyperpol/masdnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/masdnp.m`
- Signature: `dnp=masdnp(spin_system,parameters)`
- Total lines: 204

## Purpose

Magic angle spinning DNP simulation, returning the rotor period averaged steady state magnetization. This function takes a lot of inspiration from the code donated by Fred Mentink, please ci- te Fred's papers if you are using it. Syntax: dnp=masdnp(spin_system,parameters)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 56-57: Convert microwave frequency into offset; implemented by `parameters.offset=parameters.mw_frq-spin_system.inter.magnet*spin(parameters.spins{1})/(2*pi)`.
- Lines 59-61: Microwave operator; implemented by `Hmw=(operator(spin_system,'L+',parameters.spins{1})+ operator(spin_system,'L-',parameters.spins{1}))/2`.
- Lines 63-64: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 66-67: Thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 69-71: Get the averaging grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid '.mat'])`.
- Lines 73-75: Shut up and inform the user; implemented by `report(spin_system,['powder average being computed over ' num2str(numel(sph_grid.weights)) ' orientations.'])`.
- Lines 81-82: Get the answer going; implemented by `dnp=0`.
- Lines 84-85: Loop over the grid weights; implemented by `parfor n=1:numel(sph_grid.weights)`.
- Lines 87-88: Set the current orientation; implemented by `localpar=parameters`.
- Lines 94-95: Rotor stack generation; implemented by `H=rotor_stack(spin_system,localpar,'esr')`.
- Lines 98-99: Rotor period integration; implemented by `nsteps=numel(H); P=speye(size(H{1}))`.
- Lines 104-105: Effective rotor period Liouvillian; implemented by `L_eff=1i*localpar.rate*logm(P)`.
- Lines 107-108: Evolve for the equilibration time; implemented by `rho_st=evolution(spin_system,L_eff,[],rho_eq,parameters.mw_time,1,'final')`.
- Lines 110-111: Rotor period trajectory; implemented by `rho=zeros([numel(rho_eq) nsteps],'like',1i); rho(:,1)=rho_st`.
- Lines 116-117: Enhancement factor; implemented by `Hz_dnp=mean(localpar.coil'*rho)`.

### Control flow inferred from the code

- Line 76: conditional branch on `~parameters.verbose`.
- Line 85: `parfor` loop over `n=1:numel(sph_grid.weights)`.
- Line 100: `for` loop over `k=1:nsteps`.
- Line 112: `for` loop over `k=2:nsteps`.

### Key state/data transformations

- Lines 57: computes `parameters.offset` using `parameters.offset=parameters.mw_frq-spin_system.inter.magnet*spin(parameters.spins{1})/(2*pi)`.
- Lines 60-61: computes `Hmw` using `Hmw=(operator(spin_system,'L+',parameters.spins{1})+ operator(spin_system,'L-',parameters.spins{1}))/2`.
- Lines 64: computes `R` using `R=relaxation(spin_system)`.
- Lines 67: computes `rho_eq` using `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 70-71: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid '.mat'])`.
- Lines 78: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 82: computes `dnp` using `dnp=0`.
- Lines 88: computes `localpar` using `localpar=parameters`.
- Lines 89-91: computes `localpar.orientation` using `localpar.orientation=[sph_grid.alphas(n) sph_grid.betas(n) sph_grid.gammas(n)]`.
- Lines 92: computes `localpar.masframe` using `localpar.masframe='rotor'; localpar.rframes={}`.
- Lines 95: computes `H` using `H=rotor_stack(spin_system,localpar,'esr')`.
- Lines 96: computes `nsteps` using `nsteps=numel(H); dt=1/(nsteps*localpar.rate)`.
- Lines 101: computes `P` using `P=propagator(spin_system,H{k}+localpar.mw_pwr*Hmw+1i*R,dt)*P`.
- Lines 105: computes `L_eff` using `L_eff=1i*localpar.rate*logm(P)`.
- Lines 108: computes `rho_st` using `rho_st=evolution(spin_system,L_eff,[],rho_eq,parameters.mw_time,1,'final')`.
- Lines 111: computes `rho` using `rho=zeros([numel(rho_eq) nsteps],'like',1i); rho(:,1)=rho_st`.
- Lines 113: computes `rho(:,k)` using `rho(:,k)=step(spin_system,H{k}+localpar.mw_pwr*Hmw+1i*R,rho(:,k-1),dt)`.
- Lines 117: computes `Hz_dnp` using `Hz_dnp=mean(localpar.coil'*rho)`.

### Local helper functions

- Line 126: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})`.
  - Representative operation: `error('this function is only available in Liouville space.')`.

## Parameters / inputs

- parameters.spins -the spins to microwave
- parameters.rate -spinning rate, Hz
- parameters.axis -spinning axis direction vector.
- parameters.max_rank -rotor discretization grid rank
- parameters.mw_pwr -microwave power, rad/s
- parameters.mw_frq -microwave frequency, Hz
- parameters.mw_time -microwave irradiation duration
- before the average magnetistion
- is computed, seconds
- parameters.grid -the name of the spherical avera-
- ging grid
- parameters.coil -detection state
- parameters.verbose -set this to 1 to enable diag-
- nostic output

## Outputs

- dnp -enhancement of the user-specified state relative to
- the thermal equilibrium
- Note: increase the rotor rank and the spherical grid size until
- the answer stops changing. You will likely need huge values
- for both parameters.
- Note: this function must be called directly, without a context
- wrapper.

## Implementation structure

- Magic angle spinning DNP simulation, returning the rotor period
- averaged steady state magnetization. This function takes a lot
- of inspiration from the code donated by Fred Mentink, please ci-
- te Fred's papers if you are using it. Syntax:
- dnp=masdnp(spin_system,parameters)
- parameters.spins - the spins to microwave
- parameters.rate - spinning rate, Hz
- parameters.axis - spinning axis direction vector.
- parameters.max_rank - rotor discretization grid rank
- parameters.mw_pwr - microwave power, rad/s
- parameters.mw_frq - microwave frequency, Hz
- parameters.mw_time - microwave irradiation duration

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `operator()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `load()`, `report()`, `num2str()`, `rotor_stack()`, `speye()`, `propagator()`, `logm()`, `evolution()`, `rho()`.
