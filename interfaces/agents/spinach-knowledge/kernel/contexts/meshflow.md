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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 60-61: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters)`.
- Lines 63-64: Get diffusion and flow generator; implemented by `F=flow_gen(spin_system,parameters)`.
- Lines 66-67: Get problem dimensions; implemented by `spc_dim=spin_system.mesh.vor.ncells`.
- Lines 74-76: Build Hamiltonian superoperator out of phantoms; implemented by `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 84-86: Build relaxation superoperator out of phantoms; implemented by `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 94-96: Build kinetics superoperator out of phantoms; implemented by `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 104-105: Spin-independent spatial motion; implemented by `F=polyadic({{F,opium(spn_dim,1)}})`.
- Lines 107-108: Dummy gradients for now; implemented by `G=polyadic({{0}})`.
- Lines 110-111: Inflate polyadic objects if necessary; implemented by `if ~ismember('polyadic',spin_system.sys.enable)`.
- Lines 116-117: Build the initial condition out of phantoms; implemented by `parameters.rho0=0`.
- Lines 123-124: Build the detection state out of phantoms; implemented by `parameters.coil=0`.
- Lines 130-131: Call the pulse sequence; implemented by `answer=pulse_sequence(spin_system,parameters,H,R,K,G,F)`.

### Control flow inferred from the code

- Line 77: `for` loop over `n=1:numel(parameters.H_ph)`.
- Line 87: `for` loop over `n=1:numel(parameters.R_ph)`.
- Line 97: `for` loop over `n=1:numel(parameters.K_ph)`.
- Line 111: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.
- Line 118: `for` loop over `n=1:numel(parameters.rho0_ph)`.
- Line 125: `for` loop over `n=1:numel(parameters.coil_ph)`.

### Key state/data transformations

- Lines 64: computes `F` using `F=flow_gen(spin_system,parameters)`.
- Lines 67: computes `spc_dim` using `spc_dim=spin_system.mesh.vor.ncells`.
- Lines 68: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1); problem_dim=spc_dim*spn_dim`.
- Lines 72: computes `parameters.spc_dim` using `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 75-76: computes `polyadic_cores` using `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 78: computes `polyadic_cores{n}` using `polyadic_cores{n}=cell(1,2)`.
- Lines 79: computes `polyadic_cores{n}{1}` using `polyadic_cores{n}{1}=spdiags(parameters.H_ph{n}(:),0,spc_dim,spc_dim)`.
- Lines 80: computes `polyadic_cores{n}{2}` using `polyadic_cores{n}{2}=parameters.H_op{n}`.
- Lines 82: computes `H` using `H=polyadic(polyadic_cores)`.
- Lines 92: computes `R` using `R=polyadic(polyadic_cores)`.
- Lines 102: computes `K` using `K=polyadic(polyadic_cores)`.
- Lines 108: computes `G` using `G=polyadic({{0}})`.
- Lines 117: computes `parameters.rho0` using `parameters.rho0=0`.
- Lines 124: computes `parameters.coil` using `parameters.coil=0`.
- Lines 131: computes `answer` using `answer=pulse_sequence(spin_system,parameters,H,R,K,G,F)`.

### Local helper functions

- Line 136: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters)`.
  - Representative operation: `if ~isfield(spin_system,'mesh')`.
  - Representative operation: `error('mesh information is missing from the spin_system structure.')`.

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
