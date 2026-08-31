# kernel/contexts/imaging.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/imaging.m`
- Signature: `answer=imaging(spin_system,pulse_sequence,parameters)`
- Total lines: 281

## Purpose

Fokker-Planck imaging simulation context. Generates the Hamiltonian, the relaxation superoperator, the kinetics superoperator, the Fokker- Planck spatial dynamics generator (including diffusion and flow), gra- dient operators, and passes all of that to the pulse sequence, which should be supplied as a handle. Syntax: answer=imaging(spin_system,pulse_sequence,parameters)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 86-87: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 89-90: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 92-93: Set NMR assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 95-96: Call Spinach to build Hamiltonian superoperator; implemented by `H=hamiltonian(spin_system)`.
- Lines 98-99: Process channel offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 101-102: Call Spinach to build kinetics superoperator; implemented by `K=kinetics(spin_system)`.
- Lines 104-105: Get problem dimensions; implemented by `spc_dim=prod(parameters.npts); spn_dim=size(H,1); problem_dim=spc_dim*spn_dim`.
- Lines 111-113: Build relaxation superoperator out of phantoms; implemented by `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 121-122: Decide the initial state; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 124-125: Build out of phantoms; implemented by `parameters.rho0=0`.
- Lines 134-135: Use initial state as received; implemented by `report(spin_system,'initial state supplied by the user.')`.
- Lines 139-140: Decide the coil state; implemented by `if ~isfield(parameters,'coil')`.
- Lines 142-143: Build out of phantoms; implemented by `parameters.coil=0`.
- Lines 152-153: Use coil state as received; implemented by `report(spin_system,'coil state supplied by the user.')`.
- Lines 157-158: Put the same spin dynamics and chemistry into every voxel; implemented by `H=polyadic({{opium(spc_dim,1),H}})`.
- Lines 161-162: Get gradient operators; implemented by `G=g2fplanck(spin_system,parameters)`.
- Lines 164-165: Get diffusion and hydrodynamics generators; implemented by `F=v2fplanck(spin_system,parameters)`.
- Lines 167-168: Inflate polyadic objects; implemented by `if ~ismember('polyadic',spin_system.sys.enable)`.

### Control flow inferred from the code

- Line 114: `for` loop over `n=1:numel(parameters.rlx_ph)`.
- Line 122: conditional branch on `~isfield(parameters,'rho0')`.
- Line 126: `for` loop over `n=1:numel(parameters.rho0_ph)`.
- Line 140: conditional branch on `~isfield(parameters,'coil')`.
- Line 144: `for` loop over `n=1:numel(parameters.coil_ph)`.
- Line 168: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.

### Key state/data transformations

- Lines 93: computes `spin_system` using `spin_system=assume(spin_system,'nmr')`.
- Lines 96: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 102: computes `K` using `K=kinetics(spin_system)`.
- Lines 105: computes `spc_dim` using `spc_dim=prod(parameters.npts); spn_dim=size(H,1); problem_dim=spc_dim*spn_dim`.
- Lines 109: computes `parameters.spc_dim` using `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 112-113: computes `polyadic_cores` using `polyadic_cores={{sparse([],[],[],spc_dim,spc_dim), sparse([],[],[],spn_dim,spn_dim)}}`.
- Lines 115: computes `polyadic_cores{n}` using `polyadic_cores{n}=cell(1,2)`.
- Lines 116: computes `polyadic_cores{n}{1}` using `polyadic_cores{n}{1}=spdiags(parameters.rlx_ph{n}(:),0,spc_dim,spc_dim)`.
- Lines 117: computes `polyadic_cores{n}{2}` using `polyadic_cores{n}{2}=parameters.rlx_op{n}`.
- Lines 119: computes `R` using `R=polyadic(polyadic_cores)`.
- Lines 125: computes `parameters.rho0` using `parameters.rho0=0`.
- Lines 143: computes `parameters.coil` using `parameters.coil=0`.
- Lines 162: computes `G` using `G=g2fplanck(spin_system,parameters)`.
- Lines 165: computes `F` using `F=v2fplanck(spin_system,parameters)`.
- Lines 175: computes `answer` using `answer=pulse_sequence(spin_system,parameters,H,R,K,G,F)`.

### Local helper functions

- Line 180: `grumble()` — `function grumble(spin_system,parameters)`. Enforce relaxation phantom specification
  - Representative operation: `if (~isfield(parameters,'rlx_ph'))||(~isfield(parameters,'rlx_op'))`.
  - Representative operation: `error('relaxation phantom must be provided in parameters.rlx_ph and parameters.rlx_op')`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.u -X components of the velocity vectors
- for each point in the sample, m/s
- parameters.v -Y components of the velocity vectors
- for each point in the sample, m/s
- parameters.w -Z components of the velocity vectors
- for each point in the sample, m/s
- parameters.diff -diffusion coefficient or 3x3 tensor, m^2/s
- for situations when this parameter is the
- same in every voxel
- parameters.dxx -Cartesian components of the diffusion
- parameters.dxy tensor for each voxel of the sample
- ...
- parameters.dzz
- parameters.dims -dimensions of the 3D box, meters
- parameters.npts -number of points in each dimension
- of the 3D box
- parameters.deriv -{'fourier'} uses Fourier diffe-
- rentiation matrices; {'period',n}
- requests n-point central finite-
- difference matrices with periodic
- boundary conditions
- Three types of phantoms must be specified. The relaxation theory phantom
- contains relaxation superoperators and their coefficients in each voxel,
- specified in the following way:
- parameters.rlx_ph={Ph1,Ph2,...,PhN}
- parameters.rlx_op={R1,R2,...,RN}
- where PhN have the same dimension as the sample voxel grid and RN are re-
- laxation superoperators. The initial condition phantom reflects the fact
- that different voxels might start off in a different spin state. It must
- be specified in the following way:
- parameters.rho0_ph={Ph1,Ph2,...,PhN}
- parameters.rho0_st={rho1,rho2,...,rhoN}
- where PhN have the same dimension as the sample voxel grid and rhoN are
- spin states obtained from state() function. The detection state phantom
- reflects the fact that different voxels might be detected at different
- angles and with different sensitivity. It must be specified in the follo-
- wing way:
- parameters.coil_ph={Ph1,Ph2,...,PhN}
- parameters.coil_st={rho1,rho2,...,rhoN}
- where PhN have the same dimension as the sample voxel grid and rhoN are
- spin states obtained from state() function.

## Outputs

- This function returns whatever the pulse sequence returns.
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].

## Implementation structure

- Fokker-Planck imaging simulation context. Generates the Hamiltonian,
- the relaxation superoperator, the kinetics superoperator, the Fokker-
- Planck spatial dynamics generator (including diffusion and flow), gra-
- dient operators, and passes all of that to the pulse sequence, which
- should be supplied as a handle. Syntax:
- answer=imaging(spin_system,pulse_sequence,parameters)
- pulse_sequence - pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.u -X components of the velocity vectors
- for each point in the sample, m/s
- parameters.v -Y components of the velocity vectors

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `grumble()`, `assume()`, `hamiltonian()`, `frqoffset()`, `kinetics()`, `report()`, `num2str()`, `spdiags()`, `polyadic()`, `isfield()`, `opium()`, `g2fplanck()`, `v2fplanck()`, `ismember()`, `complex()`.
