# experiments/cp_contact_hard.m

- Signature: `contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)`

## Purpose

Cross-polarisation experiment in the rotating frame. Applies an ideal pi/2 pulse using the specified operators, then evolves the system with the specified spin-lock terms added to the Liovilli- an. The contact curve is returned. Syntax: contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.irr_powers -a matrix containing the values
- of the spin-lock nutation fre-
- quency on each channel (rows)
- at each time slice (cols), Hz
- parameters.irr_opers -a cell array of spin operators
- corresponding to the spin-lock
- on each channel
- parameters.exc_opers -a cell array of spin operators
- for the ideal pi/2 excitation
- pulse (same flip angle on all
- channels)
- parameters.time_steps -a vector of time slice durati-
- ons, seconds
- parameters.rho0 -initial state vector
- parameters.coil -detection state vector
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- contact_curve -contact curve detected on the coil
- state specified in parameters.coil

## Implementation structure

- Cross-polarisation experiment in the rotating frame. Applies an
- ideal pi/2 pulse using the specified operators, then evolves the
- system with the specified spin-lock terms added to the Liovilli-
- an. The contact curve is returned. Syntax:
- contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)
- parameters.irr_powers -a matrix containing the values
- of the spin-lock nutation fre-
- quency on each channel (rows)
- at each time slice (cols), Hz
- parameters.irr_opers -a cell array of spin operators
- corresponding to the spin-lock
- on each channel
