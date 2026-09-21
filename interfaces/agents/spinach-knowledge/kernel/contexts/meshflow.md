# kernel/contexts/meshflow.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/meshflow.m`
- Signature: `answer=meshflow(spin_system,pulse_sequence,parameters)`
- Total lines: 173

## Purpose

First draft of the magnetohydrodynamics context for microfluidic simu- lations. Generates evolution generators and passes them on to the pul- se sequence function, which should be supplied as a handle. Syntax: answer=meshflow(spin_system,pulse_sequence,parameters)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- The following phantoms must be specified: hamiltonian, relaxation, ki-
- netics, initial condition, detection state. Operator phantoms must be
- specified in the following way:
- parameters.R_ph={Ph1,Ph2,...,PhN}
- parameters.R_op={R1,R2,...,RN}
- where PhN have the same dimension as the sample voxel grid and RN are
- relaxation superoperators. Likewise for the following:
- parameters.K_ph, parameters.K_op
- parameters.H_ph, parameters.H_op
- The initial condition phantom reflects the fact that different voxels
- might start off in a different spin state. It must be specified in the
- following way:
- parameters.rho0_ph={Ph1,Ph2,...,PhN}
- parameters.rho0_st={rho1,rho2,...,rhoN}
- where PhN have the same dimension as the sample voxel grid and rhoN are
- spin states obtained from state() function.
- The detection state phantom reflects the fact that different voxels mi-
- ght be detected at different angles and with different sensitivity. It
- must be specified in the following way:
- parameters.coil_ph={Ph1,Ph2,...,PhN}
- parameters.coil_st={rho1,rho2,...,rhoN}
- where PhN have the same dimension as the sample voxel grid and rhoN are
- spin states obtained from state() function.
- parameters.* -additional subfields may be required by your
- pulse sequence -check its documentation page

## Outputs

- This function returns whatever the pulse sequence returns.

## Implementation structure

- First draft of the magnetohydrodynamics context for microfluidic simu-
- lations. Generates evolution generators and passes them on to the pul-
- se sequence function, which should be supplied as a handle. Syntax:
- answer=meshflow(spin_system,pulse_sequence,parameters)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- The following phantoms must be specified: hamiltonian, relaxation, ki-
- netics, initial condition, detection state. Operator phantoms must be
- specified in the following way:
- parameters.R_ph={Ph1,Ph2,...,PhN}
- parameters.R_op={R1,R2,...,RN}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `flow_gen()`, `report()`, `num2str()`, `spdiags()`, `polyadic()`, `opium()`, `ismember()`, `inflate()`, `pulse_sequence()`, `isfield()`.
