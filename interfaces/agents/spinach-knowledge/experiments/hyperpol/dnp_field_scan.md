# experiments/hyperpol/dnp_field_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/dnp_field_scan.m`
- Signature: `dnp=dnp_field_scan(spin_system,parameters,H,R,K)`
- Total lines: 185

## Purpose

Magnetic field scan steady-state DNP experiment. Returns the steady-state population of the user-specified state as a fun- ction of magnetic field. Syntax: dnp=dnp_field_scan(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 61-62: Damp the trace direction and check; implemented by `switch spin_system.bas.formalism`.
- Lines 66-67: Unit state carries the trace; implemented by `R(1,1)=-mean(abs(diag(R)))`.
- Lines 71-72: Trace direction is the stretched unit matrix; implemented by `dim=sqrt(size(R,1)); u0=speye(dim); u0=u0(:)`.
- Lines 83-84: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 86-87: Preallocate the answer; implemented by `dnp=zeros([numel(parameters.fields) size(parameters.coil,2)],'like',1i)`.
- Lines 89-90: Add microwave terms to the Liouvillian; implemented by `L=L+2*pi*parameters.mw_pwr*parameters.mw_oper`.
- Lines 92-93: Add offset terms to the Liouvillian; implemented by `L=L-2*pi*parameters.mw_frq*parameters.ez_oper`.
- Lines 95-96: Precompute the right hand side; implemented by `b=R*parameters.rho0`.
- Lines 98-99: Loop over the fields; implemented by `parfor n=1:numel(parameters.fields)`.
- Lines 101-102: Calculate electron offset frequency; implemented by `omega=-parameters.fields(n)*spin('E')`.
- Lines 104-105: Add electron offset terms to the Liouvillian; implemented by `L_current=L+omega*parameters.ez_oper`.
- Lines 107-108: Get the steady state DNP; implemented by `if strcmp(parameters.method,'gmres')`.
- Lines 110-111: Make sure L is sparse;; implemented by `L_current=sparse(L_current)`.
- Lines 113-114: Update the preconditioner; implemented by `[M1,M2]=ilu(-1i*L_current,struct('type','crout','droptol',1e-6))`.
- Lines 116-117: Run using preconditioned GMRES; implemented by `dnp(n,:)=parameters.coil'*gmres(-1i*L_current,b,10,1e-10,numel(parameters.rho0),M1,M2)`.
- Lines 121-122: Run using Matlab's backslash; implemented by `dnp(n,:)=parameters.coil'*(-1i*L_current\b)`.
- Lines 126-127: Complain and bomb out; implemented by `error('unknown solver.')`.

### Control flow inferred from the code

- Line 62: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`.
- Line 77: conditional branch on `Rc>1e9, error('R must be non-singular.'); end`.
- Line 78: conditional branch on `any(ismember({'redfield','naka-zwan'},spin_system.rlx.theories))&&`.
- Line 99: `parfor` loop over `n=1:numel(parameters.fields)`.
- Line 108: conditional branch on `strcmp(parameters.method,'gmres')`.

### Key state/data transformations

- Lines 67: computes `R(1,1)` using `R(1,1)=-mean(abs(diag(R)))`.
- Lines 72: computes `dim` using `dim=sqrt(size(R,1)); u0=speye(dim); u0=u0(:)`.
- Lines 73: computes `R` using `R=R-(mean(abs(diag(R)))/dim)*(u0*u0')`.
- Lines 76: computes `Rc` using `Rc=condest(R)`.
- Lines 84: computes `L` using `L=H+1i*R+1i*K`.
- Lines 87: computes `dnp` using `dnp=zeros([numel(parameters.fields) size(parameters.coil,2)],'like',1i)`.
- Lines 96: computes `b` using `b=R*parameters.rho0`.
- Lines 102: computes `omega` using `omega=-parameters.fields(n)*spin('E')`.
- Lines 105: computes `L_current` using `L_current=L+omega*parameters.ez_oper`.
- Lines 114: computes `[M1,M2]` using `[M1,M2]=ilu(-1i*L_current,struct('type','crout','droptol',1e-6))`.
- Lines 117: computes `dnp(n,:)` using `dnp(n,:)=parameters.coil'*gmres(-1i*L_current,b,10,1e-10,numel(parameters.rho0),M1,M2)`.

### Local helper functions

- Line 136: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.mw_pwr -microwave power, Hz
- parameters.mw_frq -microwave frequency offset from
- the free electron frequency at
- the reference B0 field, Hz
- parameters.fields -a vector of magnetic field off-
- sets from the reference B0 field,
- Tesla
- parameters.rho0 -equilibrium state at the reference
- B0 field
- parameters.coil -coil state vector or a horizon-
- tal stack thereof
- parameters.mw_oper -microwave irradiation operator
- parameters.ez_oper -Lz operator on the electrons
- parameters.method -'backslash' to use Matlab's
- linear equation solver, 'gmres'
- to use ILU preconditioned GMRES
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- dnp -an array of steady state expectation values for
- the states specified in parameters.coil at each
- of the fields supplied
- Note: the relaxation superoperator should NOT be thermalized
- for this type of calculation.
- Note: thermal equilibrium state and relaxation superoperator are
- assumed to be the same at all fields in the sweep -DO NOT
- USE with broad magnetic field sweep experiments.

## Implementation structure

- Magnetic field scan steady-state DNP experiment. Returns the
- steady-state population of the user-specified state as a fun-
- ction of magnetic field. Syntax:
- dnp=dnp_field_scan(spin_system,parameters,H,R,K)
- parameters.mw_pwr - microwave power, Hz
- parameters.mw_frq - microwave frequency offset from
- the free electron frequency at
- the reference B0 field, Hz
- parameters.fields - a vector of magnetic field off-
- sets from the reference B0 field,
- Tesla
- parameters.rho0 - equilibrium state at the reference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `condest()`, `any()`, `ismember()`, `spin()`, `strcmp()`, `ilu()`, `dnp()`, `gmres()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `isvector()`, `ischar()`.
