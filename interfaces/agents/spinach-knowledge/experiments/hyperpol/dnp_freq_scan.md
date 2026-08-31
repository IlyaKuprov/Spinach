# experiments/hyperpol/dnp_freq_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/dnp_freq_scan.m`
- Signature: `dnp=dnp_freq_scan(spin_system,parameters,H,R,K)`
- Total lines: 282

## Purpose

Microwave frequency scan steady-state DNP experiment. Returns the steady-state population of the user-specified states as a function of microwave irradiation frequency. Syntax: dnp=dnp_freq_scan(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 61-62: Damp the trace direction and check; implemented by `switch spin_system.bas.formalism`.
- Lines 66-67: Unit state carries the trace; implemented by `R(1,1)=-mean(abs(diag(R)))`.
- Lines 71-72: Trace direction is the stretched unit matrix; implemented by `dim=sqrt(size(R,1)); u0=speye(dim); u0=u0(:)`.
- Lines 83-84: Preallocate the answer; implemented by `dnp=zeros([numel(parameters.mw_frq) size(parameters.coil,2)],'like',1i)`.
- Lines 88-89: Fokker-Planck path; implemented by `case {'fp-backs','fp-gmres'}`.
- Lines 91-92: Translate frequency offsets; implemented by `carrier_freq=-spin('E')*(parameters.g_ref/spin_system.tols.freeg)*spin_system.inter.magnet`.
- Lines 95-96: Get problem dimensions; implemented by `spc_dim=parameters.nphases; spn_dim=size(H,1)`.
- Lines 98-99: Get phases and Fourier derivative operator; implemented by `[phases,d_dphi]=fourdif(spc_dim,1)`.
- Lines 101-102: Build the phase turning operator; implemented by `PT=kron(d_dphi,speye(spn_dim,spn_dim))`.
- Lines 104-105: Replicate static Hamiltonian; implemented by `H=kron(speye(spc_dim,spc_dim),H)`.
- Lines 107-108: Replicate relaxation and kinetics; implemented by `R=kron(speye(spc_dim,spc_dim),R)`.
- Lines 111-112: Build microwave operator; implemented by `cosines=spdiags(cos(phases),0,spc_dim,spc_dim)`.
- Lines 115-116: Precompute the right hand side; implemented by `b=R*kron(ones(spc_dim,1),parameters.rho0)`.
- Lines 118-119: MDCS diagnostics; implemented by `parallel_profiler_start`.
- Lines 121-122: Loop in parallel over the frequencies; implemented by `parfor n=1:numel(parameters.mw_frq)`.
- Lines 124-125: Build the evolution generator; implemented by `L=sparse(H+1i*R+1i*K+MW+1i*parameters.mw_frq(n)*PT)`.
- Lines 127-128: Get the steady state DNP; implemented by `switch parameters.method`.

### Control flow inferred from the code

- Line 62: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `{'fp-backs','fp-gmres'}`, `'fp-gmres'`, `'fp-backs'`.
- Line 77: conditional branch on `Rc>1e9, error('R must be non-singular.'); end`.
- Line 78: conditional branch on `any(ismember({'redfield','naka-zwan'},spin_system.rlx.theories))&&`.
- Line 86: dispatches on `parameters.method`; cases `{'fp-backs','fp-gmres'}`, `'fp-gmres'`, `'fp-backs'`.
- Line 122: `parfor` loop over `n=1:numel(parameters.mw_frq)`.
- Line 128: dispatches on `parameters.method`; cases `'fp-gmres'`, `'fp-backs'`.
- Line 172: `parfor` loop over `n=1:numel(parameters.mw_frq)`.
- Line 178: dispatches on `parameters.method`; cases `'lvn-gmres'`, `'lvn-backs'`.

### Key state/data transformations

- Lines 67: computes `R(1,1)` using `R(1,1)=-mean(abs(diag(R)))`.
- Lines 72: computes `dim` using `dim=sqrt(size(R,1)); u0=speye(dim); u0=u0(:)`.
- Lines 73: computes `R` using `R=R-(mean(abs(diag(R)))/dim)*(u0*u0')`.
- Lines 76: computes `Rc` using `Rc=condest(R)`.
- Lines 84: computes `dnp` using `dnp=zeros([numel(parameters.mw_frq) size(parameters.coil,2)],'like',1i)`.
- Lines 92: computes `carrier_freq` using `carrier_freq=-spin('E')*(parameters.g_ref/spin_system.tols.freeg)*spin_system.inter.magnet`.
- Lines 93: computes `parameters.mw_frq` using `parameters.mw_frq=carrier_freq+parameters.mw_frq`.
- Lines 96: computes `spc_dim` using `spc_dim=parameters.nphases; spn_dim=size(H,1)`.
- Lines 99: computes `[phases,d_dphi]` using `[phases,d_dphi]=fourdif(spc_dim,1)`.
- Lines 102: computes `PT` using `PT=kron(d_dphi,speye(spn_dim,spn_dim))`.
- Lines 105: computes `H` using `H=kron(speye(spc_dim,spc_dim),H)`.
- Lines 109: computes `K` using `K=kron(speye(spc_dim,spc_dim),K)`.
- Lines 112: computes `cosines` using `cosines=spdiags(cos(phases),0,spc_dim,spc_dim)`.
- Lines 113: computes `MW` using `MW=parameters.mw_pwr*kron(cosines,parameters.mw_oper)`.
- Lines 116: computes `b` using `b=R*kron(ones(spc_dim,1),parameters.rho0)`.
- Lines 125: computes `L` using `L=sparse(H+1i*R+1i*K+MW+1i*parameters.mw_frq(n)*PT)`.
- Lines 133: computes `big_inv` using `big_inv=gmres(-1i*L,b,10,1e-10,numel(b),L)`.
- Lines 148: computes `dnp(n,:)` using `dnp(n,:)=parameters.coil'*mean(reshape(big_inv,[spn_dim spc_dim]),2)`.

### Local helper functions

- Line 215: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.mw_pwr -microwave power, rad/s
- parameters.mw_frq -row vector of microwave frequ-
- ency offsets (rad/s) relative
- to the reference g-factor
- parameters.g_ref -reference g-factor around which
- frequency offsets are specified
- parameters.rho0 -thermal equilibrium state
- parameters.coil -coil state vector or a horizon-
- tal stack thereof
- parameters.mw_oper -microwave irradiation operator
- parameters.ez_oper -Lz operator on the electrons
- parameters.method -calculation method: 'fp-backs',
- 'fp-gmres', 'lvn-backs', or
- 'lvn-gmres'
- parameters.nphases -number of microwave phase grid
- points for the Fokker-Planck path
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- dnp -an array of steady state expectation values for
- the states specified in parameters.coil at each
- of the microwave frequencies supplied
- Note: the relaxation superoperator must NOT be thermalized for
- this type of calculation (inter.equilibrium='zero').

## Implementation structure

- Microwave frequency scan steady-state DNP experiment. Returns the
- steady-state population of the user-specified states as a function
- of microwave irradiation frequency. Syntax:
- dnp=dnp_freq_scan(spin_system,parameters,H,R,K)
- parameters.mw_pwr - microwave power, rad/s
- parameters.mw_frq - row vector of microwave frequ-
- ency offsets (rad/s) relative
- to the reference g-factor
- parameters.g_ref - reference g-factor around which
- frequency offsets are specified
- parameters.rho0 - thermal equilibrium state
- parameters.coil - coil state vector or a horizon-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `condest()`, `any()`, `ismember()`, `spin()`, `fourdif()`, `spdiags()`, `gmres()`, `dnp()`, `ilu()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `isscalar()`.
