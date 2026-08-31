# experiments/overtone/overtone_cp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/overtone/overtone_cp.m`
- Signature: `spectrum=overtone_cp(spin_system,parameters,H,R,K)`
- Total lines: 176

## Purpose

Cross-polarization overtone experiment. Syntax: spectrum=overtone_cp(spin_system,parameters,H,R,K)

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

- Lines 62-63: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 65-66: Get the overtone frequency; implemented by `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 68-69: Project pulse operators; implemented by `Nx=kron(speye(parameters.spc_dim),parameters.Nx)`.
- Lines 72-73: Choose the method; implemented by `switch parameters.method`.
- Lines 77-78: Get the frequency; implemented by `omega=2*pi*ovt_frq-2*pi*parameters.rf_frq`.
- Lines 80-82: Get the pulse Hamiltonian; implemented by `pulseop=average(spin_system,parameters.rf_pwr(1)*Nx/2,H+parameters.rf_pwr(2)*Hx, parameters.rf_pwr(1)*Nx/2,omega,'matrix_log')`.
- Lines 84-85: Apply the pulse; implemented by `parameters.rho0=propagator(spin_system,pulseop,parameters.rf_dur)*parameters.rho0`.
- Lines 89-90: Get the frequency; implemented by `omega=ovt_frq-parameters.rf_frq`.
- Lines 92-95: Apply the pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,H+parameters.rf_pwr(2)*Hx,Nx,0*Nx, parameters.rho0,omega,parameters.rf_pwr(1), parameters.rf_dur,-pi/2,2,'expm')`.
- Lines 99-100: Call the acquisition; implemented by `spectrum=overtone_a(spin_system,parameters,H,R,K)`.

### Control flow inferred from the code

- Line 73: dispatches on `parameters.method`; cases `'average'`, `'fplanck'`.

### Key state/data transformations

- Lines 66: computes `ovt_frq` using `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 69: computes `Nx` using `Nx=kron(speye(parameters.spc_dim),parameters.Nx)`.
- Lines 70: computes `Hx` using `Hx=kron(speye(parameters.spc_dim),parameters.Hx)`.
- Lines 78: computes `omega` using `omega=2*pi*ovt_frq-2*pi*parameters.rf_frq`.
- Lines 81-82: computes `pulseop` using `pulseop=average(spin_system,parameters.rf_pwr(1)*Nx/2,H+parameters.rf_pwr(2)*Hx, parameters.rf_pwr(1)*Nx/2,omega,'matrix_log')`.
- Lines 85: computes `parameters.rho0` using `parameters.rho0=propagator(spin_system,pulseop,parameters.rf_dur)*parameters.rho0`.
- Lines 100: computes `spectrum` using `spectrum=overtone_a(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 105: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.spins overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim Fokker-Planck spatial dimension
- parameters.method pulse simulation method, either
- 'average' or 'fplanck'
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Nx X Zeeman operator on the
- quadrupolar nucleus
- parameters.Hx X Zeeman operator on the
- spin-1/2 nucleus
- parameters.rf_frq spin-lock frequency offset from
- the overtone frequency on the
- quadrupolar nucleus, Hz
- parameters.rf_pwr a vector of spin-lock powers on
- the quadrupolar nucleus (first
- element) and the spin-1/2 nucleus
- (second element), rad/s
- parameters.rf_dur spin-lock pulse duration, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the resulting spectrum
- Notes: relaxation must be present in the system dynamics, or the
- matrix inversion in overtone_a function call would fail to
- converge. The relaxation matrix must *not* be thermalised.

## Implementation structure

- Cross-polarization overtone experiment. Syntax:
- spectrum=overtone_cp(spin_system,parameters,H,R,K)
- parameters.spins overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim Fokker-Planck spatial dimension
- parameters.method pulse simulation method, either
- 'average' or 'fplanck'
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `speye()`, `average()`, `propagator()`, `shaped_pulse_af()`, `overtone_a()`, `ismatrix()`, `isfield()`, `elseif()`, `iscell()`, `ismember()`.
