# experiments/esr_dipolar/eseem.m

- Signature: `fid=eseem(spin_system,parameters,H,R,K)`

## Purpose

ESEEM pulse sequence with ideal hard pulses. Syntax: fid=eseem(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.npoints number of points to be computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -time domaing signal whose Fourier transform is the
- ESEEM spectrum

## Implementation structure

- ESEEM pulse sequence with ideal hard pulses. Syntax:
- fid=eseem(spin_system,parameters,H,R,K)
- parameters.npoints number of points to be computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
