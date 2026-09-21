# experiments/overtone/overtone_pa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/overtone/overtone_pa.m`
- Signature: `spectrum=overtone_pa(spin_system,parameters,H,R,K)`
- Total lines: 180

## Purpose

Overtone soft pulse-acquire experiment. Syntax: spectrum=overtone_pa(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Overtone experiment implementations. These routines excite or detect high-order quadrupolar transitions and therefore combine non-secular quadrupolar terms, MAS or field effects, and specialised detection pathways.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.spc_dim Fokker-Planck spatial dimension
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Lx X Zeeman operator on the
- quadrupolar nucleus
- parameters.rf_frq pulse frequency offset from
- the overtone frequency on the
- quadrupolar nucleus, Hz
- parameters.rf_pwr pulse power on the quadrupolar
- nucleus, rad/s
- parameters.rf_dur pulse duration, seconds
- parameters.method 'average' uses the average Hamil-
- tonian theory, 'fplanck' uses
- Fokker-Planck formalism
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the matrix
- inversion operation in the overtone_a call would fail. The rela-
- xation superoperator R must *not* be thermalised.

## Implementation structure

- Overtone soft pulse-acquire experiment. Syntax:
- spectrum=overtone_pa(spin_system,parameters,H,R,K)
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz around the overtone frequency
- parameters.npoints number of points in the spectrum
- parameters.spins overtone-active nucleus, specified as a
- single-element cell array
- parameters.spc_dim Fokker-Planck spatial dimension
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.Lx X Zeeman operator on the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `speye()`, `average()`, `propagator()`, `shaped_pulse_af()`, `overtone_a()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ismember()`.
