# experiments/hyperpol/solid_effect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/solid_effect.m`
- Signature: `answer=solid_effect(spin_system,parameters)`
- Total lines: 219

## Purpose

Solid effect DNP experiment, computed using the large-scale formalism described in (http://dx.doi.org/10.1039/C2CP23233B). The system is re- stricted to one electron and one nucleus type, but the number of nuc- lei may be very large. Syntax: answer=solid_effect(spin_system,parameters)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `find()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 60-61: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 63-64: Build H+ and add microwave term; implemented by `report(spin_system,'solid_effect: building H+ ')`.
- Lines 68-69: Build H0; implemented by `report(spin_system,'solid_effect: building H0 ')`.
- Lines 73-74: Build H-and add microwave term; implemented by `report(spin_system,'solid_effect: building H- ')`.
- Lines 78-79: Decide how to proceed; implemented by `switch parameters.theory`.
- Lines 83-85: Build the exact Hamiltonian; implemented by `H=Hp+H0+Hm+parameters.nuclear_frq*operator(spin_system,'Lz','nuclei')+ -parameters.nuclear_frq*operator(spin_system,'Lz','E')`.
- Lines 89-90: Build an average Hamiltonian; implemented by `H=average(spin_system,Hp,H0,Hm,parameters.nuclear_frq,parameters.theory)`.
- Lines 94-95: Set state option; implemented by `if isfield(parameters,'coil')`.
- Lines 97-98: Use user-specified states; implemented by `coils=parameters.coil`.
- Lines 102-103: Set detection states to Lz on every spin; implemented by `coils=cell(1,spin_system.comp.nspins)`.
- Lines 111-112: Set the initial state to thermal equilibrium; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 115-116: Decide the output; implemented by `switch parameters.calc_type`.
- Lines 120-121: Get thermalized relaxation superoperator; implemented by `spin_system.rlx.equilibrium='IME'`.
- Lines 124-126: Run the simulation; implemented by `answer=evolution(spin_system,H+1i*R,coils,rho_eq,parameters.time_step, parameters.n_steps,'multichannel')`.
- Lines 130-131: Get unthermalized relaxation superoperator; implemented by `spin_system.rlx.equilibrium='zero'; R=relaxation(spin_system,[0 0 0])`.
- Lines 133-134: Remove unit state singularity; implemented by `rho_unit=unit_state(spin_system); R=R-kron(rho_unit,rho_unit')`.
- Lines 136-137: Run the simulation; implemented by `answer=evolution(spin_system,H+1i*R,coils,-R*rho_eq,[],[],'total')`.
- Lines 141-142: Get thermalized relaxation superoperator; implemented by `spin_system.rlx.equilibrium='IME'; R=relaxation(spin_system)`.

### Control flow inferred from the code

- Line 79: dispatches on `parameters.theory`; cases `'exact'`.
- Line 95: conditional branch on `isfield(parameters,'coil')`.
- Line 104: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 116: dispatches on `parameters.calc_type`; cases `'time_dependence'`, `'steady_state'`, `'trajectory'`.

### Key state/data transformations

- Lines 65: computes `[Hp,Q]` using `[Hp,Q]=hamiltonian(assume(spin_system,'se_dnp_h+'))`.
- Lines 66: computes `Hp` using `Hp=Hp+orientation(Q,[0 0 0])+0.25*parameters.mw_pwr*operator(spin_system,'L-','E')`.
- Lines 70: computes `[H0,Q]` using `[H0,Q]=hamiltonian(assume(spin_system,'se_dnp_h0'))`.
- Lines 71: computes `H0` using `H0=H0+orientation(Q,[0 0 0])`.
- Lines 75: computes `[Hm,Q]` using `[Hm,Q]=hamiltonian(assume(spin_system,'se_dnp_h-'))`.
- Lines 76: computes `Hm` using `Hm=Hm+orientation(Q,[0 0 0])+0.25*parameters.mw_pwr*operator(spin_system,'L+','E')`.
- Lines 84-85: computes `H` using `H=Hp+H0+Hm+parameters.nuclear_frq*operator(spin_system,'Lz','nuclei')+ -parameters.nuclear_frq*operator(spin_system,'Lz','E')`.
- Lines 98: computes `coils` using `coils=parameters.coil`.
- Lines 105: computes `coils{n}` using `coils{n}=state(spin_system,{'Lz'},{n})`.
- Lines 112: computes `[I,Q]` using `[I,Q]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 113: computes `rho_eq` using `rho_eq=equilibrium(spin_system,I,Q,[0 0 0])`.
- Lines 121: computes `spin_system.rlx.equilibrium` using `spin_system.rlx.equilibrium='IME'`.
- Lines 122: computes `R` using `R=relaxation(spin_system,[0 0 0])`.
- Lines 125-126: computes `answer` using `answer=evolution(spin_system,H+1i*R,coils,rho_eq,parameters.time_step, parameters.n_steps,'multichannel')`.
- Lines 134: computes `rho_unit` using `rho_unit=unit_state(spin_system); R=R-kron(rho_unit,rho_unit')`.

### Local helper functions

- Line 153: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.mw_pwr -microwave power in rad/s
- parameters.theory -level of theory. Set to 'exact' for
- the electron rotating frame calcula-
- tion or to any of the following six
- options for the average Hamiltonian
- theory calculation on top of the el-
- ectron + nuclear rotating frame:
- 'ah_first_order', 'ah_second_order',
- 'ah_third_order', 'kb_first_order',
- 'kb_second_order', 'kb_third_order'.
- See average.m function for the mea-
- ning of these options.
- parameters.nuclear_frq -nuclear Zeeman frequency in rad/s
- parameters.calc_type -set to 'time_dependence' to get the
- time dependence of the longitudinal
- magnetization and to 'steady_state'
- to get the asymptotic longitudinal
- magnetization.
- parameters.time_step -if 'time_dependence' is set in the
- calc_type parameter, sets the time
- step, seconds.
- parameters.n_steps -if 'time_dependence' is set in the
- calc_type parameter, sets the num-
- ber of time steps.

## Outputs

- answer -with the 'time_dependence' calculation type, the
- function returns the observables detected using
- the coil states specified at each point in time;
- with the 'steady_state' option specified, the
- function returns the steady state values detec-
- ted using the coil states specified.
- Note: this function generates its own Liouvillian and should be
- called directly, without a context wrapper.

## Implementation structure

- Solid effect DNP experiment, computed using the large-scale formalism
- described in (http://dx.doi.org/10.1039/C2CP23233B). The system is re-
- stricted to one electron and one nucleus type, but the number of nuc-
- lei may be very large. Syntax:
- answer=solid_effect(spin_system,parameters)
- parameters.mw_pwr -microwave power in rad/s
- parameters.theory -level of theory. Set to 'exact' for
- the electron rotating frame calcula-
- tion or to any of the following six
- options for the average Hamiltonian
- theory calculation on top of the el-
- ectron + nuclear rotating frame:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `hamiltonian()`, `assume()`, `orientation()`, `operator()`, `average()`, `isfield()`, `state()`, `cell2mat()`, `equilibrium()`, `relaxation()`, `evolution()`, `unit_state()`, `ismember()`, `cellfun()`.
