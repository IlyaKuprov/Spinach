# experiments/esr_dipolar/deer_3p_hard_deer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_3p_hard_deer.m`
- Signature: `deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)`
- Total lines: 204

## Purpose

Three-pulse DEER pulse sequence. Idealized hard pulses are used, each pulse only affects its specific electron or transition, de- pending on the pulse operators supplied. Syntax: deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 70-71: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 73-74: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 76-77: Compute detailed output if required; implemented by `if strcmp(parameters.output,'detailed')`.
- Lines 79-80: Apply a hard 90-degree pulse; implemented by `rho_hard=step(spin_system,parameters.ex_hard,parameters.rho0,pi/2)`.
- Lines 82-84: Return the free induction decay; implemented by `deer.hard_pulse_fid=evolution(spin_system,L,parameters.coil_pump+parameters.coil_prob,rho_hard, 1/parameters.spectrum_sweep,parameters.spectrum_nsteps,'observable')`.
- Lines 86-87: Apply a selective 90-degree pulse on the probe spin; implemented by `rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2)`.
- Lines 89-91: Return the free induction decay on the probe spin; implemented by `deer.prob_pulse_fid=evolution(spin_system,L,parameters.coil_prob,rho,1/parameters.spectrum_sweep, parameters.spectrum_nsteps,'observable')`.
- Lines 93-94: Apply a selective 90-degree pulse on the pump spin; implemented by `rho=step(spin_system,parameters.ex_pump,parameters.rho0,pi/2)`.
- Lines 96-98: Return the free induction decay on the pump spin; implemented by `deer.pump_pulse_fid=evolution(spin_system,L,parameters.coil_pump,rho,1/parameters.spectrum_sweep, parameters.spectrum_nsteps,'observable')`.
- Lines 102-103: First pulse; implemented by `rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2)`.
- Lines 105-106: Evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,parameters.stepsize,parameters.nsteps,'trajectory')`.
- Lines 108-109: Second pulse; implemented by `rho_stack=step(spin_system,parameters.ex_pump,rho_stack,pi)`.
- Lines 111-112: Evolution; implemented by `rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.stepsize,parameters.nsteps,'refocus')`.
- Lines 114-115: Third pulse; implemented by `rho_stack=step(spin_system,parameters.ex_prob,rho_stack,pi)`.
- Lines 117-118: Evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.stepsize,parameters.nsteps,'final')`.
- Lines 120-121: Observation; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-hilb')`.

### Control flow inferred from the code

- Line 77: conditional branch on `strcmp(parameters.output,'detailed')`.
- Line 121: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.

### Key state/data transformations

- Lines 74: computes `L` using `L=H+1i*R+1i*K`.
- Lines 80: computes `rho_hard` using `rho_hard=step(spin_system,parameters.ex_hard,parameters.rho0,pi/2)`.
- Lines 83-84: computes `deer.hard_pulse_fid` using `deer.hard_pulse_fid=evolution(spin_system,L,parameters.coil_pump+parameters.coil_prob,rho_hard, 1/parameters.spectrum_sweep,parameters.spectrum_nsteps,'observable')`.
- Lines 87: computes `rho` using `rho=step(spin_system,parameters.ex_prob,parameters.rho0,pi/2)`.
- Lines 90-91: computes `deer.prob_pulse_fid` using `deer.prob_pulse_fid=evolution(spin_system,L,parameters.coil_prob,rho,1/parameters.spectrum_sweep, parameters.spectrum_nsteps,'observable')`.
- Lines 97-98: computes `deer.pump_pulse_fid` using `deer.pump_pulse_fid=evolution(spin_system,L,parameters.coil_pump,rho,1/parameters.spectrum_sweep, parameters.spectrum_nsteps,'observable')`.
- Lines 106: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,parameters.stepsize,parameters.nsteps,'trajectory')`.
- Lines 112: computes `rho_stack(:,end:-1:1)` using `rho_stack(:,end:-1:1)=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.stepsize,parameters.nsteps,'refocus')`.
- Lines 122: computes `deer.deer_trace` using `deer.deer_trace=cellfun(@(x)(trace(parameters.coil_prob'*x)),rho_stack)/norm(parameters.coil_prob,'fro')`.

### Local helper functions

- Line 130: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.rho0 initial state
- parameters.coil_prob detection state on probe spin
- parameters.stepsize increment time for the pump pulse
- sandwich
- parameters.nsteps number of steps for the pump pulse
- sandwich
- parameters.ex_prob excitation operators to be used for
- parameters.ex_pump the probe and pump electron respec-
- tively.
- parameters.output 'brief' returns just the DEER trace,
- 'detailed' also returns excitation
- profiles and the EPR spectrum.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- If 'detailed' is selected as the output option, the following pa-
- rameters are also required:
- parameters.ex_hard hard pulse excitation operator
- parameters.spectrum_sweep sweep width of the EPR spectrum, Hz
- parameters.spectrum_nsteps number of time steps in the FID
- parameters.coil_pump detection state on pump spin

## Outputs

- deer.hard_pulse_fid -('detailed') free induction decay
- after a non-selective ideal pulse
- deer.prob_pulse_fid -('detailed') free induction decay
- after just the the probe pulse
- deer.pump_pulse_fid -('detailed') free induction decay
- after just the pump pulse
- deer.deer_trace -DEER signal
- Note: hard pulses are only appropriate for spin-1/2 systems; for
- higher spin systems transition selective pulse operators
- must be supplied.

## Implementation structure

- Three-pulse DEER pulse sequence. Idealized hard pulses are used,
- each pulse only affects its specific electron or transition, de-
- pending on the pulse operators supplied. Syntax:
- deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)
- parameters.rho0 initial state
- parameters.coil_prob detection state on probe spin
- parameters.stepsize increment time for the pump pulse
- sandwich
- parameters.nsteps number of steps for the pump pulse
- parameters.ex_prob excitation operators to be used for
- parameters.ex_pump the probe and pump electron respec-
- tively.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `step()`, `evolution()`, `rho_stack()`, `cellfun()`, `ismatrix()`, `all()`, `isfield()`, `ischar()`, `ismember()`.
