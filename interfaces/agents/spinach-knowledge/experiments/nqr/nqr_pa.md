# experiments/nqr/nqr_pa.m

- Signature: `spectrum=nqr_pa(spin_system,parameters,H,R,K)`

## Purpose

Nuclear quadrupole resonance soft pulse-acquire experiment. Idealised ac- quisition with infinite bandwidth is done. Syntax: spectrum=nqr_pa(spin_system,parameters,H,R,K)

## Physical / mathematical content

- NQR experiment implementations. These pulse sequences work in quadrupolar-dominated regimes with little or no Zeeman interaction and focus on nutation, free evolution, and transition detection in electric-field-gradient frames.

## Numerical / algorithmic content

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
