# kernel/utilities/tolerances.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/tolerances.m`
- Signature: `[spin_system,sys]=tolerances(spin_system,sys)`
- Total lines: 519

## Purpose

Tolerances and fundamental constants. Sets various accuracy cut-offs, constants and tolerances used by Spinach kernel. Syntax: spin_system=tolerances(spin_system,sys)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(spin_system,sys)`.
- Lines 35-36: Interaction tensor clean-up tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')`.
- Lines 54-55: Liouvillian matrix element zero tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')`.
- Lines 73-74: Relaxation superoperator element zero tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'rlx_zero')`.
- Lines 92-93: Propagator matrix element zero tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')`.
- Lines 111-112: Subspace population tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'subs_drop')`.
- Lines 130-131: Irrep population tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')`.
- Lines 149-150: Steady state convergence tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'stst_tol')`.
- Lines 168-169: ZTE sample length; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')`.
- Lines 187-188: ZTE zero track tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')`.
- Lines 206-207: ZTE state vector density threshold; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')`.
- Lines 221-222: Proximity tolerance for dipolar couplings (TODO: replace with energy tolerance); implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')`.
- Lines 240-241: Krylov method switchover; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'krylov_tol')`.
- Lines 251-252: Basis printing hush tolerance; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')`.
- Lines 262-263: Subspace bundle size; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'merge_dim')`.
- Lines 273-274: Sparse algebra tolerance on density; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'dense_matrix')`.
- Lines 284-285: Sparse algebra tolerance on dimension; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'small_matrix')`.
- Lines 295-296: Relative accuracy of the elements of Redfield superoperator; implemented by `if isfield(sys,'tols')&&isfield(sys.tols,'rlx_integration')`.

### Control flow inferred from the code

- Line 36: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')`.
- Line 55: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')`.
- Line 74: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'rlx_zero')`.
- Line 93: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')`.
- Line 112: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'subs_drop')`.
- Line 131: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')`.
- Line 150: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'stst_tol')`.
- Line 169: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')`.
- Line 188: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')`.
- Line 207: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')`.
- Line 222: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')`.
- Line 241: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'krylov_tol')`.
- Line 252: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')`.
- Line 263: conditional branch on `isfield(sys,'tols')&&isfield(sys.tols,'merge_dim')`.

### Key state/data transformations

- Lines 37: computes `spin_system.tols.inter_cutoff` using `spin_system.tols.inter_cutoff=sys.tols.inter_cutoff; sys.tols=rmfield(sys.tols,'inter_cutoff')`.
- Lines 56: computes `spin_system.tols.liouv_zero` using `spin_system.tols.liouv_zero=sys.tols.liouv_zero; sys.tols=rmfield(sys.tols,'liouv_zero')`.
- Lines 75: computes `spin_system.tols.rlx_zero` using `spin_system.tols.rlx_zero=sys.tols.rlx_zero; sys.tols=rmfield(sys.tols,'rlx_zero')`.
- Lines 94: computes `spin_system.tols.prop_chop` using `spin_system.tols.prop_chop=sys.tols.prop_chop; sys.tols=rmfield(sys.tols,'prop_chop')`.
- Lines 113: computes `spin_system.tols.subs_drop` using `spin_system.tols.subs_drop=sys.tols.subs_drop; sys.tols=rmfield(sys.tols,'subs_drop')`.
- Lines 132: computes `spin_system.tols.irrep_drop` using `spin_system.tols.irrep_drop=sys.tols.irrep_drop; sys.tols=rmfield(sys.tols,'irrep_drop')`.
- Lines 151: computes `spin_system.tols.stst_tol` using `spin_system.tols.stst_tol=sys.tols.stst_tol; sys.tols=rmfield(sys.tols,'stst_tol')`.
- Lines 170: computes `spin_system.tols.zte_nsteps` using `spin_system.tols.zte_nsteps=sys.tols.zte_nsteps; sys.tols=rmfield(sys.tols,'zte_nsteps')`.
- Lines 189: computes `spin_system.tols.zte_tol` using `spin_system.tols.zte_tol=sys.tols.zte_tol; sys.tols=rmfield(sys.tols,'zte_tol')`.
- Lines 208: computes `spin_system.tols.zte_maxden` using `spin_system.tols.zte_maxden=sys.tols.zte_maxden; sys.tols=rmfield(sys.tols,'zte_maxden')`.
- Lines 223: computes `spin_system.tols.prox_cutoff` using `spin_system.tols.prox_cutoff=sys.tols.prox_cutoff; sys.tols=rmfield(sys.tols,'prox_cutoff')`.
- Lines 242: computes `spin_system.tols.krylov_tol` using `spin_system.tols.krylov_tol=sys.tols.krylov_tol; sys.tols=rmfield(sys.tols,'krylov_tol')`.
- Lines 253: computes `spin_system.tols.basis_hush` using `spin_system.tols.basis_hush=sys.tols.basis_hush; sys.tols=rmfield(sys.tols,'basis_hush')`.
- Lines 264: computes `spin_system.tols.merge_dim` using `spin_system.tols.merge_dim=sys.tols.merge_dim; sys.tols=rmfield(sys.tols,'merge_dim')`.
- Lines 275: computes `spin_system.tols.dense_matrix` using `spin_system.tols.dense_matrix=sys.tols.dense_matrix; sys.tols=rmfield(sys.tols,'dense_matrix')`.
- Lines 286: computes `spin_system.tols.small_matrix` using `spin_system.tols.small_matrix=sys.tols.small_matrix; sys.tols=rmfield(sys.tols,'small_matrix')`.
- Lines 297: computes `spin_system.tols.rlx_integration` using `spin_system.tols.rlx_integration=sys.tols.rlx_integration; sys.tols=rmfield(sys.tols,'rlx_integration')`.
- Lines 316: computes `spin_system.tols.dP_method` using `spin_system.tols.dP_method=sys.tols.dP_method; sys.tols=rmfield(sys.tols,'dP_method')`.

### Local helper functions

- Line 390: `grumble()` — `function grumble(spin_system,sys)`.
  - Representative operation: `if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')`.
  - Representative operation: `if (~isnumeric(sys.tols.inter_cutoff))||(~isreal(sys.tols.inter_cutoff))|| (~isscalar(sys.tols.inter_cutoff))||(sys.tols.inter_cutoff<0)`.

## Parameters / inputs

- spin_system -Spinach system description object
- sys -system specification object described
- in the input preparation section of
- the manual

## Outputs

- spin_system -updated system description object
- sys -system specification structure with
- the tolerance substructure parsed out
- Notes: direct calls and modifications to this function are discouraged:
- the accuracy settings should be modified by setting the sys.tols
- structure, see the input preparation manual.

## Implementation structure

- Tolerances and fundamental constants. Sets various accuracy cut-offs,
- constants and tolerances used by Spinach kernel. Syntax:
- spin_system=tolerances(spin_system,sys)
- spin_system - Spinach system description object
- sys - system specification object described
- in the input preparation section of
- the manual
- spin_system - updated system description object
- sys - system specification structure with
- the tolerance substructure parsed out
- the accuracy settings should be modified by setting the sys.tols
- structure, see the input preparation manual.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rmfield()`, `report()`, `pad()`, `below()`, `num2str()`, `ismember()`, `eps()`, `than()`, `inf()`, `constant()`, `setdiff()`, `fieldnames()`, `isscalar()`, `ischar()`.
