# kernel/optimcon/optimcon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/optimcon.m`
- Signature: `spin_system=optimcon(spin_system,control)`
- Total lines: 1480

## Purpose

Validates optimal control options and updates the spin system object. Syntax: spin_system=optimcon(spin_system,control)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `check_hermiticity()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Consistency check; implemented by `grumble(spin_system,control)`.
- Lines 46-47: Delete the previous control structure; implemented by `if isfield(spin_system,'control')`.
- Lines 52-53: Show the banner; implemented by `banner(spin_system,'optimcon')`.
- Lines 55-56: Process fidelity type; implemented by `if isfield(control,'fidelity')`.
- Lines 58-59: Input validation; implemented by `if ~ischar(control.fidelity)`.
- Lines 63-64: Absorb the input; implemented by `spin_system.control.fidelity=control.fidelity`.
- Lines 69-70: Default is Re(<targ|P|init>); implemented by `spin_system.control.fidelity='real'`.
- Lines 74-75: Inform the user; implemented by `switch spin_system.control.fidelity`.
- Lines 79-81: Real part of the overlap; implemented by `report(spin_system,[pad('Fidelity measure, range [-1,+1]',60) pad('Re(<target|rho(T)>)',20)])`.
- Lines 85-87: Imaginary part of the overlap; implemented by `report(spin_system,[pad('Fidelity measure, range [-1,+1]',60) pad('Im(<target|rho(T)>)',20)])`.
- Lines 91-93: Absolute square of the overlap; implemented by `report(spin_system,[pad('Fidelity measure, range [0,+1]',60) pad('|<target|rho(T)>|^2',20)])`.
- Lines 97-98: Complain and bomb out; implemented by `error('control.fidelity can be ''real'', ''imag'', or ''square''.')`.
- Lines 102-103: Process integrator type; implemented by `if isfield(control,'integrator')`.
- Lines 105-106: Input validation; implemented by `if ~ischar(control.integrator)`.
- Lines 113-114: Absorb integrator type; implemented by `spin_system.control.integrator=control.integrator`.
- Lines 119-120: Default is piecewise-constant; implemented by `spin_system.control.integrator='rectangle'`.
- Lines 124-126: Inform the user; implemented by `report(spin_system,[pad('Equation of motion integrator',60) spin_system.control.integrator])`.
- Lines 128-129: Process optimisation method; implemented by `if isfield(control,'method')`.

### Control flow inferred from the code

- Line 47: conditional branch on `isfield(spin_system,'control')`.
- Line 56: conditional branch on `isfield(control,'fidelity')`.
- Line 59: conditional branch on `~ischar(control.fidelity)`.
- Line 75: dispatches on `spin_system.control.fidelity`; cases `'real'`, `'imag'`, `'square'`.
- Line 103: conditional branch on `isfield(control,'integrator')`.
- Line 106: conditional branch on `~ischar(control.integrator)`.
- Line 109: conditional branch on `~ismember(control.integrator,{'rectangle','trapezium'})`.
- Line 129: conditional branch on `isfield(control,'method')`.
- Line 132: conditional branch on `(~ischar(control.method))||(~ismember(control.method,{'lbfgs','rbfgs',`.
- Line 142: conditional branch on `strcmp(spin_system.control.method,'newton')&&`.
- Line 148: conditional branch on `strcmp(spin_system.control.method,'goodwin')&&`.
- Line 169: conditional branch on `isfield(control,'isotopes')`.
- Line 172: conditional branch on `(~iscell(control.isotopes))||(~all(cellfun(@ischar,control.isotopes(:))))`.
- Line 175: `for` loop over `n=1:numel(control.isotopes)`.

### Key state/data transformations

- Lines 48: computes `spin_system` using `spin_system=rmfield(spin_system,'control')`.
- Lines 64: computes `spin_system.control.fidelity` using `spin_system.control.fidelity=control.fidelity`.
- Lines 65: computes `control` using `control=rmfield(control,'fidelity')`.
- Lines 114: computes `spin_system.control.integrator` using `spin_system.control.integrator=control.integrator`.
- Lines 138: computes `spin_system.control.method` using `spin_system.control.method=control.method`.
- Lines 164-166: computes `herm_gens` using `herm_gens=ismember(spin_system.bas.formalism, {'zeeman-hilb','zeeman-wavef'})|| strcmp(spin_system.control.method,'goodwin')`.
- Lines 182: computes `spin_system.control.isotopes` using `spin_system.control.isotopes=control.isotopes`.
- Lines 208: computes `spin_system.control.channels` using `spin_system.control.channels=control.channels(:)`.
- Lines 232: computes `spin_system.control.ncontrols` using `spin_system.control.ncontrols=numel(control.operators)`.
- Lines 241-242: computes `spin_system.control.operators{n}` using `spin_system.control.operators{n}=clean_up(spin_system,control.operators{n}, spin_system.tols.liouv_zero)`.
- Lines 254-255: computes `spin_system.control.cc_comm` using `spin_system.control.cc_comm=cell(spin_system.control.ncontrols, spin_system.control.ncontrols)`.
- Lines 256-257: computes `spin_system.control.cc_comm_idx` using `spin_system.control.cc_comm_idx=false(spin_system.control.ncontrols, spin_system.control.ncontrols)`.
- Lines 264-265: computes `spin_system.control.cc_comm{n,m}` using `spin_system.control.cc_comm{n,m}=comm(spin_system.control.operators{n}, spin_system.control.operators{m})`.
- Lines 268: computes `comm_norm` using `comm_norm=norm(spin_system.control.cc_comm{n,m},1)`.
- Lines 269: computes `spin_system.control.cc_comm_idx(n,m)` using `spin_system.control.cc_comm_idx(n,m)=(comm_norm<spin_system.tols.liouv_zero)`.
- Lines 275: computes `n_commute` using `n_commute=(nnz(spin_system.control.cc_comm_idx)-numel(spin_system.control.operators))/2`.
- Lines 328: computes `spin_system.control.rho_init` using `spin_system.control.rho_init=control.rho_init`.
- Lines 388: computes `spin_system.control.rho_targ` using `spin_system.control.rho_targ=control.rho_targ`.

### Local helper functions

- Line 1444: `check_hermiticity()` — `function check_hermiticity(generators,generator_kind)`. Over the supplied generators
  - Representative operation: `for n=1:numel(generators)`.
  - Representative operation: `norm_a=cheap_norm(generators{n}-generators{n}')`.
- Line 1465: `grumble()` — `function grumble(spin_system,control)`. Whenever anyone accuses some person of being 'unfeeling' he means that that person is just. He means that that person has no causeless emotions
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a Spinach data structure.')`.

## Parameters / inputs

- spin_system -primary Spinach data structure,
- created by create.m and updated
- by basis.m functions
- control -control data structure described
- in detail in the online manual

## Outputs

- spin_system -updated Spinach data structure
- Note: this function freezes the optimisation problem. The ensemble
- case catalog is built here, and the complete frozen problem
- is published to the parallel pool workers exactly once, as a
- parallel.pool.Constant held in spin_system.control.invari-
- ants. Heavy invariants -the drift generators, the control
- operators, the offset operators, the control commutators,
- and the Bloch-Siegert response operators -are then removed
- from the returned structure, and their names are recorded
- in spin_system.control.frozen_fields. All other control
- fields stay live: ensemble() re-sends them to the workers
- at every evaluation, and they may be overwritten between
- optimisations. Changes to the ensemble composition, the
- operators, or the generators require a fresh optimcon()
- call.

## Implementation structure

- Validates optimal control options and updates the spin system
- object. Syntax:
- spin_system=optimcon(spin_system,control)
- spin_system -primary Spinach data structure,
- created by create.m and updated
- by basis.m functions
- control -control data structure described
- in detail in the online manual
- spin_system -updated Spinach data structure
- Note: this function freezes the optimisation problem. The ensemble
- case catalog is built here, and the complete frozen problem
- is published to the parallel pool workers exactly once, as a

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rmfield()`, `report()`, `banner()`, `ischar()`, `pad()`, `rho()`, `ismember()`, `strcmp()`, `iscell()`, `all()`, `cellfun()`, `isvector()`, `any()`, `check_hermiticity()`.
