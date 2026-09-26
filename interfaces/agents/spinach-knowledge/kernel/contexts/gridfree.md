# kernel/contexts/gridfree.m

- Signature: `answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil- lian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. Syntax: answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Parameters / inputs

- pulse_sequence -a function handle to one of the pulse sequences
- located in the experiments directory
- assumptions -is a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters -a structure with the following subfields:
- .rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.
- .axis -spinning axis, given as a normalized
- 3-element vector
- .spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- .offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- .max_rank -maximum D-function rank to retain in
- the solution (increase till conver-
- gence is achieved, approximately
- equal to the number of spinning si-
- debands in the spectrum)
- .tau_c -correlation times (in seconds) for rotational
- diffusion. Single number for isotropic rotati-
- onal diffusion, and a symmetric positive defi-
- nite 3x3 correlation time tensor for anisotro-
- pic rotational diffusion; the rotational dif-
- fusion tensor is inv(6*tau_c).
- .* -additional subfields may be required by your
- pulse sequence -check its documentation page
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- this context function returns the powder average of whatever it
- is that the pulse sequence returns
- Note: the choice of the Wigner D function rank truncation level depends on
- on the spinning rate (the slower the spinning, the greater ranks are
- required).
- Note: rotational correlation times for SLE go into parameters.tau_c, not
- inter.tau_c (the latter is only used by the Redfield theory module).
- Note: the state projector assumes a powder --single crystal MAS is not
- currently supported.
- Note: perturbative corrections to the rotating frame transformation are
- not supported -use singlerot.m if you need them.

## Implementation structure

- Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil-
- lian superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle. Syntax:
- answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -a function handle to one of the pulse sequences
- located in the experiments directory
- assumptions -is a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters -a structure with the following subfields:
- .rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.
