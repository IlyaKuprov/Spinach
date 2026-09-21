# experiments/spin_chem/rydmr_exp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spin_chem/rydmr_exp.m`
- Signature: `answer=rydmr_exp(spin_system,parameters,H,R,K)`
- Total lines: 179

## Purpose

Singlet-singlet RYDMR experiment with exponential recombination function (http://dx.doi.org/10.1080/00268979809483134). Syntax: A=rydmr_exp(spin_system,parameters,H,R,K) where H is the Hamiltonian commutation superoperator in zero ex- ternal field, R is the relaxation superoperator and K is the che- mical kinetics superoperator.

## Physical / mathematical content

- Spin-chemistry experiment implementations. These routines couple spin evolution to chemical kinetics, radical-pair recombination, exchange, and spin-selective reaction channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.fields -row vector of field values, Tesla; the
- primary magnet field should be set to
- sys.magnet=1 for normalisation purposes
- parameters.rates -row vector of singlet recombination
- rate constants, Hz
- parameters.electrons -numbers identifying the two electrons
- in the isotope list, e.g. [1 2]
- parameters.needs -must contain 'zeeman_op', this is an
- instruction to the kernel to provide a
- separate Zeeman operator for field sweep
- purposes

## Outputs

- A -a matrix of singlet yields with dimensions
- matching the sizes of parameters.rates and
- parameters.fields
- Note: exponential recombination kinetics is built into this func-
- tion, do not combine with inter.chem.rp_rates parameter.

## Implementation structure

- Singlet-singlet RYDMR experiment with exponential recombination
- function (http://dx.doi.org/10.1080/00268979809483134). Syntax:
- A=rydmr_exp(spin_system,parameters,H,R,K)
- where H is the Hamiltonian commutation superoperator in zero ex-
- ternal field, R is the relaxation superoperator and K is the che-
- mical kinetics superoperator.
- parameters.fields - row vector of field values, Tesla; the
- primary magnet field should be set to
- sys.magnet=1 for normalisation purposes
- parameters.rates - row vector of singlet recombination
- rate constants, Hz
- parameters.electrons -numbers identifying the two electrons

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `singlet()`, `ind2sub()`, `fields()`, `rates()`, `speye()`, `answer()`, `evolution()`, `unit_state()`, `hdot()`, `expmint()`, `ismatrix()`, `all()`, `specification()`, `isfield()`, `iscell()`.
