# experiments/hyperpol/dnp_freq_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/dnp_freq_scan.m`
- Signature: `dnp=dnp_freq_scan(spin_system,parameters,H,R,K)`
- Total lines: 282

## Purpose

Microwave frequency scan steady-state DNP experiment. Returns the steady-state population of the user-specified states as a function of microwave irradiation frequency. Syntax: dnp=dnp_freq_scan(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.mw_pwr -microwave power, rad/s
- parameters.mw_frq -row vector of microwave frequ-
- ency offsets (rad/s) relative
- to the reference g-factor
- parameters.g_ref -reference g-factor around which
- frequency offsets are specified
- parameters.rho0 -thermal equilibrium state
- parameters.coil -coil state vector or a horizon-
- tal stack thereof
- parameters.mw_oper -microwave irradiation operator
- parameters.ez_oper -Lz operator on the electrons
- parameters.method -calculation method: 'fp-backs',
- 'fp-gmres', 'lvn-backs', or
- 'lvn-gmres'
- parameters.nphases -number of microwave phase grid
- points for the Fokker-Planck path
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- dnp -an array of steady state expectation values for
- the states specified in parameters.coil at each
- of the microwave frequencies supplied
- Note: the relaxation superoperator must NOT be thermalized for
- this type of calculation (inter.equilibrium='zero').

## Implementation structure

- Microwave frequency scan steady-state DNP experiment. Returns the
- steady-state population of the user-specified states as a function
- of microwave irradiation frequency. Syntax:
- dnp=dnp_freq_scan(spin_system,parameters,H,R,K)
- parameters.mw_pwr - microwave power, rad/s
- parameters.mw_frq - row vector of microwave frequ-
- ency offsets (rad/s) relative
- to the reference g-factor
- parameters.g_ref - reference g-factor around which
- frequency offsets are specified
- parameters.rho0 - thermal equilibrium state
- parameters.coil - coil state vector or a horizon-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `condest()`, `any()`, `ismember()`, `spin()`, `fourdif()`, `spdiags()`, `gmres()`, `dnp()`, `ilu()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `isscalar()`.
