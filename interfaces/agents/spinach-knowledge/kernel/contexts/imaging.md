# kernel/contexts/imaging.m

- Signature: `answer=imaging(spin_system,pulse_sequence,parameters)`

## Purpose

Fokker-Planck imaging simulation context. Generates the Hamiltonian, the relaxation superoperator, the kinetics superoperator, the Fokker- Planck spatial dynamics generator (including diffusion and flow), gra- dient operators, and passes all of that to the pulse sequence, which should be supplied as a handle. Syntax: answer=imaging(spin_system,pulse_sequence,parameters)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
