# experiments/nqr/nqr_pa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nqr/nqr_pa.m`
- Signature: `spectrum=nqr_pa(spin_system,parameters,H,R,K)`
- Total lines: 140

## Purpose

Nuclear quadrupole resonance soft pulse-acquire experiment. Idealised ac- quisition with infinite bandwidth is done. Syntax: spectrum=nqr_pa(spin_system,parameters,H,R,K)

## Physical / mathematical content

- NQR experiment implementations. These pulse sequences work in quadrupolar-dominated regimes with little or no Zeeman interaction and focus on nutation, free evolution, and transition detection in electric-field-gradient frames.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 54-55: Project pulse operators; implemented by `Lx=kron(speye(parameters.spc_dim),parameters.Lx)`.
- Lines 58-61: Apply the soft off-resonance pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,H+1i*R+1i*K,Lx,Ly,parameters.rho0, parameters.rf_frq,parameters.rf_pwr, parameters.rf_dur,0,2,'expm')`.
- Lines 63-64: Call frequency domain acquisition; implemented by `spectrum=slowpass(spin_system,parameters,H,R,K)`.

### Key state/data transformations

- Lines 55: computes `Lx` using `Lx=kron(speye(parameters.spc_dim),parameters.Lx)`.
- Lines 56: computes `Ly` using `Ly=kron(speye(parameters.spc_dim),parameters.Ly)`.
- Lines 59-61: computes `parameters.rho0` using `parameters.rho0=shaped_pulse_af(spin_system,H+1i*R+1i*K,Lx,Ly,parameters.rho0, parameters.rf_frq,parameters.rf_pwr, parameters.rf_dur,0,2,'expm')`.
- Lines 64: computes `spectrum` using `spectrum=slowpass(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.sweep vector with two elements giving
- the spectrum window extents, Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Lx Lx and Ly operators that go into
- parameters.Ly the RF Hamiltonian
- parameters.rf_frq RF irradiation frequency, Hz
- parameters.rf_pwr the multiplier (rad/s) in front
- of [Lx*cos(ωt)+Ly*sin(ωt)] in the
- RF Hamiltonian
- parameters.rf_dur pulse duration, seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the
- matrix inversion operation would fail to converge. The re-
- laxation matrix R should *not* be thermalised.

## Implementation structure

- Nuclear quadrupole resonance soft pulse-acquire experiment. Idealised ac-
- quisition with infinite bandwidth is done. Syntax:
- spectrum=nqr_pa(spin_system,parameters,H,R,K)
- parameters.sweep vector with two elements giving
- the spectrum window extents, Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Lx Lx and Ly operators that go into
- parameters.Ly the RF Hamiltonian
- parameters.rf_frq RF irradiation frequency, Hz
- parameters.rf_pwr the multiplier (rad/s) in front

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `shaped_pulse_af()`, `slowpass()`, `ismatrix()`, `all()`, `isfield()`, `isrow()`, `isequal()`, `isscalar()`.
