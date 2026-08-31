# experiments/overtone/overtone_dante.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/overtone/overtone_dante.m`
- Signature: `spectrum=overtone_dante(spin_system,parameters,H,R,K)`
- Total lines: 189

## Purpose

Overtone DANTE experiment with frequency-domain acquisition.

## Physical / mathematical content

- Overtone experiment implementations. These routines excite or detect high-order quadrupolar transitions and therefore combine non-secular quadrupolar terms, MAS or field effects, and specialised detection pathways.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 59-60: Get the overtone frequency; implemented by `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 62-63: Project pulse operators; implemented by `Lx=kron(speye(parameters.spc_dim),parameters.Lx)`.
- Lines 65-66: Timing parameters; implemented by `rotor_period=abs(1/parameters.rate)`.
- Lines 69-70: Bomb out if the schedule makes no sense; implemented by `if (cycle_length-parameters.pulse_dur)<0`.
- Lines 74-75: Get the pulse frequency; implemented by `omega=2*pi*ovt_frq-2*pi*parameters.rf_frq`.
- Lines 77-78: Get the pulse Hamiltonian; implemented by `pulseop=parameters.pulse_amp*Lx`.
- Lines 81-82: Precompute pulse propagator; implemented by `PP=propagator(spin_system,pulseop,parameters.pulse_dur)`.
- Lines 84-85: Precompute evolution propagator; implemented by `PE=propagator(spin_system,H,cycle_length-parameters.pulse_dur)`.
- Lines 87-88: Combine the propagators; implemented by `P=clean_up(spin_system,PE*PP,spin_system.tols.prop_chop)`.
- Lines 90-92: Apply the DANTE pulse train; implemented by `parameters.rho0=multiprop(spin_system,P,parameters.rho0, parameters.n_periods*parameters.pulse_num)`.
- Lines 94-95: Call the acquisition; implemented by `spectrum=overtone_a(spin_system,parameters,H,R,K)`.

### Control flow inferred from the code

- Line 70: conditional branch on `(cycle_length-parameters.pulse_dur)<0`.

### Key state/data transformations

- Lines 60: computes `ovt_frq` using `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 63: computes `Lx` using `Lx=kron(speye(parameters.spc_dim),parameters.Lx)`.
- Lines 66: computes `rotor_period` using `rotor_period=abs(1/parameters.rate)`.
- Lines 67: computes `cycle_length` using `cycle_length=rotor_period/parameters.pulse_num`.
- Lines 75: computes `omega` using `omega=2*pi*ovt_frq-2*pi*parameters.rf_frq`.
- Lines 78: computes `pulseop` using `pulseop=parameters.pulse_amp*Lx`.
- Lines 82: computes `PP` using `PP=propagator(spin_system,pulseop,parameters.pulse_dur)`.
- Lines 85: computes `PE` using `PE=propagator(spin_system,H,cycle_length-parameters.pulse_dur)`.
- Lines 88: computes `P` using `P=clean_up(spin_system,PE*PP,spin_system.tols.prop_chop)`.
- Lines 91-92: computes `parameters.rho0` using `parameters.rho0=multiprop(spin_system,P,parameters.rho0, parameters.n_periods*parameters.pulse_num)`.
- Lines 95: computes `spectrum` using `spectrum=overtone_a(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 100: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Syntax

```matlab
spectrum=overtone_dante(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.pulse_dur -duration of the pulse, seconds
- parameters.pulse_amp -amplitude of the pulse, rad/s
- parameters.pulse_num -number of pulses within the
- rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.Lx -X Zeeman operator on the
- quadrupolar nucleus
- parameters.rf_frq -pulse frequency offset from the
- overtone frequency, Hz
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep range, Hz
- parameters.npoints -number of acquisition points
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -overtone spectrum

## Implementation structure

- Overtone DANTE experiment with frequency-domain acquisition.
- spectrum=overtone_dante(spin_system,parameters,H,R,K)
- parameters.pulse_dur -duration of the pulse, seconds
- parameters.pulse_amp -amplitude of the pulse, rad/s
- parameters.pulse_num -number of pulses within the
- rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.Lx -X Zeeman operator on the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `speye()`, `average()`, `propagator()`, `clean_up()`, `multiprop()`, `overtone_a()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`.
