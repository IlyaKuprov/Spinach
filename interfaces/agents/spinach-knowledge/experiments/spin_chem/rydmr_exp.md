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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 47-48: Isolate the Zeeman operator; implemented by `Z=parameters.hzeeman; H=H-Z`.
- Lines 50-52: Get the two-electron singlet state; implemented by `S=singlet(spin_system,parameters.electrons(1), parameters.electrons(2))`.
- Lines 54-55: Localize input arrays; implemented by `rates=parameters.rates; fields=parameters.fields`.
- Lines 57-58: Preallocate result arrays; implemented by `N=numel(rates); M=numel(fields)`.
- Lines 61-62: Decide the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 66-67: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 69-70: Normalise the singlet; implemented by `S=S/norm(S,2)`.
- Lines 72-73: Merged parfor loop; implemented by `parfor nm=1:(N*M)`.
- Lines 75-76: Extract indices; implemented by `[n,m]=ind2sub([N M],nm)`.
- Lines 78-79: Assemble Liouvillian; implemented by `L_current=L+fields(m)*Z-1i*rates(n)*speye(size(L))`.
- Lines 81-82: Compute RYDMR; implemented by `answer(nm)=rates(n)*evolution(spin_system,L_current,S,S,[],[],'total')`.
- Lines 88-89: Normalize the singlet; implemented by `S=S/norm(S,'fro')`.
- Lines 91-92: Get the unit operator; implemented by `Id=unit_state(spin_system)`.
- Lines 100-101: Assemble and tidy up Hamiltonian; implemented by `H_curr=H+fields(m)*Z`.
- Lines 104-105: Compute integration endpoint; implemented by `t_end=10/rates(n)`.
- Lines 107-108: Compute RYDMR (rates inside for roundoff reasons); implemented by `answer(nm)=hdot(S,expmint(spin_system,H_curr,S*rates(n),H_curr+1i*rates(n)*Id,t_end))`.
- Lines 114-115: Complain and bomb out; implemented by `error('unknown formalism specification.')`.

### Control flow inferred from the code

- Line 62: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv'}`, `{'zeeman-hilb'}`.
- Line 73: `parfor` loop over `nm=1:(N*M)`.
- Line 95: `parfor` loop over `nm=1:(N*M)`.

### Key state/data transformations

- Lines 48: computes `Z` using `Z=parameters.hzeeman; H=H-Z`.
- Lines 51-52: computes `S` using `S=singlet(spin_system,parameters.electrons(1), parameters.electrons(2))`.
- Lines 55: computes `rates` using `rates=parameters.rates; fields=parameters.fields`.
- Lines 58: computes `N` using `N=numel(rates); M=numel(fields)`.
- Lines 59: computes `answer` using `answer=zeros([N*M 1],'like',1i)`.
- Lines 67: computes `L` using `L=H+1i*R+1i*K`.
- Lines 76: computes `[n,m]` using `[n,m]=ind2sub([N M],nm)`.
- Lines 79: computes `L_current` using `L_current=L+fields(m)*Z-1i*rates(n)*speye(size(L))`.
- Lines 82: computes `answer(nm)` using `answer(nm)=rates(n)*evolution(spin_system,L_current,S,S,[],[],'total')`.
- Lines 92: computes `Id` using `Id=unit_state(spin_system)`.
- Lines 101: computes `H_curr` using `H_curr=H+fields(m)*Z`.
- Lines 105: computes `t_end` using `t_end=10/rates(n)`.

### Local helper functions

- Line 125: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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
