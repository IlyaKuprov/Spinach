# experiments/pulsed_field.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pulsed_field.m`
- Signature: `answer=pulsed_field(spin_system,parameters,H,R,K) %#ok<INUSD>`
- Total lines: 203

## Purpose

Magnetisation dynamics under a time-dependent magnetic field along the Z axis of the laboratory frame with spin-phonon relaxation, as measured in pulsed-field magnetometry of molecular magnets. The field profile is replaced by a staircase; on each stair the Hamil- tonian is constant, the spin-phonon dissipator is rebuilt in the eigenbasis of that Hamiltonian, and the density matrix is propa- gated in that eigenbasis 

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 92-93: Check consistency; implemented by `grumble(spin_system,parameters,H)`.
- Lines 95-96: Put the coils into a cell array; implemented by `if iscell(parameters.coil), coils=parameters.coil; else, coils={parameters.coil}; end`.
- Lines 98-99: Preallocate the output; implemented by `nrec=floor(parameters.nsteps/parameters.nout)`.
- Lines 102-103: Remove the unit field Zeeman term supplied by the context; implemented by `H=H-parameters.hzeeman; H=(H+H')/2`.
- Lines 105-106: Thermal equilibrium at zero field as the initial state; implemented by `[V,E]=eig(full(H),'vector'); pops=exp(-spin_system.tols.hbar*(E-min(E))/(spin_system.tols.kbol*spin_system.rlx.temperature))`.
- Lines 109-110: Loop over the stairs; implemented by `for n=1:parameters.nsteps`.
- Lines 112-113: Field at the stair midpoint and the Hamiltonian on the stair; implemented by `field=parameters.field_prof((n-0.5)*dt)`.
- Lines 116-117: Eigensystem of the stair Hamiltonian and the coherent half-stair phases; implemented by `[V,E]=eig(H_curr,'vector'); phases=exp(-1i*(E-E.')*dt/2)`.
- Lines 119-120: Spin-phonon coupling operator and its thermally dressed form in the eigenbasis; implemented by `XE=V'*parameters.phonon_x*V; XE=(XE+XE')/2`.
- Lines 123-124: Dissipator as matrix products in the eigenbasis; implemented by `dissip=@(rho)-pi*(XE*(RE*rho)-(RE*rho)*XE+rho*(RE'*XE)-(XE*rho)*RE')`.
- Lines 126-127: Symmetric split step in the eigenbasis; implemented by `rho_eig=phases.*(V'*rho*V); drho=dissip(rho_eig)`.
- Lines 131-132: Record the observables and report progress; implemented by `if mod(n,parameters.nout)==0`.

### Control flow inferred from the code

- Line 96: conditional branch on `iscell(parameters.coil), coils=parameters.coil; else, coils={parameters.coil}; end`.
- Line 110: `for` loop over `n=1:parameters.nsteps`.
- Line 132: conditional branch on `mod(n,parameters.nout)==0`.
- Line 134: `for` loop over `k=1:numel(coils)`.

### Key state/data transformations

- Lines 99: computes `nrec` using `nrec=floor(parameters.nsteps/parameters.nout)`.
- Lines 100: computes `answer.t` using `answer.t=zeros(nrec,1); answer.field=zeros(nrec,1); answer.obs=zeros(nrec,numel(coils))`.
- Lines 103: computes `H` using `H=H-parameters.hzeeman; H=(H+H')/2`.
- Lines 106: computes `[V,E]` using `[V,E]=eig(full(H),'vector'); pops=exp(-spin_system.tols.hbar*(E-min(E))/(spin_system.tols.kbol*spin_system.rlx.temperature))`.
- Lines 107: computes `rho` using `rho=V*diag(pops/sum(pops))*V'; dt=parameters.timestep; nrec=0`.
- Lines 113: computes `field` using `field=parameters.field_prof((n-0.5)*dt)`.
- Lines 114: computes `H_curr` using `H_curr=H+field*parameters.hzeeman; H_curr=full((H_curr+H_curr')/2)`.
- Lines 120: computes `XE` using `XE=V'*parameters.phonon_x*V; XE=(XE+XE')/2`.
- Lines 121: computes `RE` using `RE=phonon_oper(spin_system,E,XE,parameters.phonon_i0,parameters.phonon_alpha,spin_system.rlx.temperature)`.
- Lines 124: computes `dissip` using `dissip=@(rho)-pi*(XE*(RE*rho)-(RE*rho)*XE+rho*(RE'*XE)-(XE*rho)*RE')`.
- Lines 127: computes `rho_eig` using `rho_eig=phases.*(V'*rho*V); drho=dissip(rho_eig)`.
- Lines 135: computes `answer.obs(nrec,k)` using `answer.obs(nrec,k)=real(trace(coils{k}'*rho))`.

### Local helper functions

- Line 146: `grumble()` — `function grumble(spin_system,parameters,H)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
  - Representative operation: `error('this function is only available in zeeman-hilb formalism.')`.

## Parameters / inputs

- parameters.field_prof -function handle returning the field
- in Tesla at a time in seconds
- parameters.hzeeman -Zeeman operator per Tesla, rad/s/T,
- Hilbert space, supplied by the con-
- text when 'zeeman_op' is requested
- in parameters.needs
- parameters.timestep -stair width, seconds
- parameters.nsteps -number of stairs
- parameters.coil -Hilbert space observable operator or
- a cell array of them
- parameters.phonon_x -spin-phonon coupling operator, see
- rlx_phonon.m
- parameters.phonon_i0 -phonon spectral density prefactor,
- see rlx_phonon.m
- parameters.phonon_alpha -phonon spectral density exponent,
- 1 or above, see rlx_phonon.m
- parameters.nout -number of stairs between recorded
- observable values
- H -Hamiltonian received from the context function, Hilbert
- space, containing the Zeeman term at sys.magnet=1 Tesla;
- the function removes that term and adds the field of
- each stair itself
- R -relaxation superoperator received from the context
- function; ignored, the spin-phonon dissipator is built
- here at every stair
- K -kinetics superoperator received from the context
- function; ignored

## Outputs

- answer.t -column of recording times, seconds
- answer.field -column of field values at those times, Tesla
- answer.obs -matrix of observable expectation values, one
- column per coil, at the recording times
- Note: the sequence works in zeeman-hilb formalism under the crystal
- and powder contexts, which assemble the anisotropic part of
- the Hamiltonian; the liquid context drops that part, and with
- it the crystal field of a giant spin. The powder context must
- be called with parameters.sum_up=false because the answer is
- a structure; additional rotating frames (parameters.rframes)
- are not supported because the field operator is added in the
- laboratory frame. The temperature of the phonon bath is
- inter.temperature.
- Note: sys.magnet must be 1 Tesla, so that parameters.hzeeman is
- the Zeeman operator per Tesla; the Hamiltonian received from
- the context then contains the Zeeman term at 1 Tesla, which
- this function removes before adding the field on each stair.
- The initial state is the thermal equilibrium of the field-
- free Hamiltonian at the temperature of the phonon bath.
- Note: the field on each stair is evaluated at the stair midpoint.

## Implementation structure

- Magnetisation dynamics under a time-dependent magnetic field along
- the Z axis of the laboratory frame with spin-phonon relaxation, as
- measured in pulsed-field magnetometry of molecular magnets. The
- field profile is replaced by a staircase; on each stair the Hamil-
- tonian is constant, the spin-phonon dissipator is rebuilt in the
- eigenbasis of that Hamiltonian, and the density matrix is propa-
- gated in that eigenbasis by a symmetric split: exact coherent
- phases for half a stair, the dissipative step to second order in
- the dissipator times the stair width, and the phases again. The
- dissipator times the stair width must be small; the coherent part
- is treated exactly for any stair width. The dissipator is applied
- as Hilbert space matrix products (see phonon_oper.m), so the cost

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `iscell()`, `phonon_oper()`, `dissip()`, `report()`, `num2str()`, `strcmp()`, `isfield()`, `any()`, `isscalar()`, `ishermitian()`.
