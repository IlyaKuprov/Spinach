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
