# experiments/hyperpol/dnp_time_dep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/dnp_time_dep.m`
- Signature: `answer=dnp_time_dep(spin_system,parameters,H,R,K)`
- Total lines: 122

## Purpose

Time-domain spin dynamics under microwave irradiation. Syntax: answer=dnp_time_dep(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.mw_pwr -microwave power, rad/s
- parameters.mw_off -microwave frequency offset
- from free electron,rad/s
- parameters.rho0 -thermal equilibrium state
- parameters.coil -coil state vector or a hori-
- zontal stack thereof
- parameters.mw_oper -microwave irradiation operator
- parameters.ez_oper -Lz operator on the electrons
- parameters.dt -time step, seconds
- parameters.nsteps -number of time steps
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- answer -a matrix of projections of the trajectory on
- each of the coils provided at each time step
- Note: the relaxation superoperator must be thermalised for this
- type of calculation.

## Implementation structure

- Time-domain spin dynamics under microwave irradiation. Syntax:
- answer=dnp_time_dep(spin_system,parameters,H,R,K)
- parameters.mw_pwr - microwave power, rad/s
- parameters.mw_off - microwave frequency offset
- from free electron,rad/s
- parameters.rho0 - thermal equilibrium state
- parameters.coil - coil state vector or a hori-
- zontal stack thereof
- parameters.mw_oper - microwave irradiation operator
- parameters.ez_oper - Lz operator on the electrons
- parameters.dt - time step, seconds
- parameters.nsteps - number of time steps

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `evolution()`, `ismember()`, `isfield()`.
