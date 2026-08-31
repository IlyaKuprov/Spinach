# kernel/relaxation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/relaxation.m`
- Signature: `R=relaxation(spin_system,euler_angles)`
- Total lines: 888

## Purpose

Relaxation superoperator. Syntax: R=relaxation(spin_system,euler_angles)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Set defaults; implemented by `spin_system=defaults(spin_system)`.
- Lines 44-45: Check consistency; implemented by `grumble(spin_system)`.
- Lines 47-48: Get the matrix going; implemented by `R=mprealloc(spin_system,1)`.
- Lines 50-52: Detect dissipative bosonic modes; implemented by `boson_terms=isfield(spin_system.inter,'modes')&& any([spin_system.inter.modes.damp spin_system.inter.modes.dephase]>0)`.
- Lines 54-55: Warn when mode dissipation cannot be represented; implemented by `if boson_terms&&(~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'}))`.
- Lines 60-61: Do nothing if none specified; implemented by `if isempty(spin_system.rlx.theories)&&(~boson_terms)`.
- Lines 63-64: Update the user; implemented by `report(spin_system,'relaxation superoperator set to zero.')`.
- Lines 66-67: Exit the function; implemented by `return`.
- Lines 71-72: Add extended T1/T2 model terms; implemented by `if ismember('t1_t2',spin_system.rlx.theories)`.
- Lines 74-75: Inform the user; implemented by `report(spin_system,'adding extended T1/T2 model terms ')`.
- Lines 77-78: Call the extended T1/T2 model function; implemented by `if exist('euler_angles','var')`.
- Lines 80-81: Relaxation rates may depend on the orientation; implemented by `[R1,R2]=rlx_t1_t2(spin_system,euler_angles)`.
- Lines 85-86: Orientation-independent relaxation rates; implemented by `[R1,R2]=rlx_t1_t2(spin_system)`.
- Lines 90-91: Add tothe total; implemented by `R=R+R1+R2`.
- Lines 95-96: Add Redfield terms; implemented by `if ismember('redfield',spin_system.rlx.theories)`.
- Lines 98-99: Update the user and get the timer going; implemented by `report(spin_system,'Redfield theory (rotational diffusion).')`.
- Lines 102-103: Catch zero correlation times; implemented by `for n=1:numel(spin_system.rlx.tau_c)`.
- Lines 109-110: Get the rotational basis, including the non-secular terms; implemented by `report(spin_system,'computing the lab frame Hamiltonian superoperator ')`.

### Control flow inferred from the code

- Line 55: conditional branch on `boson_terms&&(~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'}))`.
- Line 61: conditional branch on `isempty(spin_system.rlx.theories)&&(~boson_terms)`.
- Line 72: conditional branch on `ismember('t1_t2',spin_system.rlx.theories)`.
- Line 78: conditional branch on `exist('euler_angles','var')`.
- Line 96: conditional branch on `ismember('redfield',spin_system.rlx.theories)`.
- Line 103: `for` loop over `n=1:numel(spin_system.rlx.tau_c)`.
- Line 104: conditional branch on `~any(spin_system.rlx.tau_c{n},'all')`.
- Line 117: conditional branch on `isworkernode||ismember('asyredf',spin_system.sys.disable)`.
- Line 127: `for` loop over `n=1:numel(spin_system.rlx.tau_c)`.
- Line 128: conditional branch on `(1/max_rate)<spin_system.rlx.tau_c{n}`.
- Line 142: conditional branch on `ismember('naka-zwan',spin_system.rlx.theories)`.
- Line 149: `for` loop over `n=1:numel(spin_system.rlx.tau_c)`.
- Line 150: conditional branch on `~any(spin_system.rlx.tau_c{n},'all')`.
- Line 156: conditional branch on `ischar(spin_system.rlx.nz_shift)`.

### Key state/data transformations

- Lines 42: computes `spin_system` using `spin_system=defaults(spin_system)`.
- Lines 48: computes `R` using `R=mprealloc(spin_system,1)`.
- Lines 51-52: computes `boson_terms` using `boson_terms=isfield(spin_system.inter,'modes')&& any([spin_system.inter.modes.damp spin_system.inter.modes.dephase]>0)`.
- Lines 81: computes `[R1,R2]` using `[R1,R2]=rlx_t1_t2(spin_system,euler_angles)`.
- Lines 100: computes `timer_redfield` using `timer_redfield=tic`.
- Lines 111: computes `[L0,Q]` using `[L0,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 114: computes `rlx_onshell` using `rlx_onshell=true; rlx_shift=0`.
- Lines 126: computes `max_rate` using `max_rate=max(abs(diag(R)))`.
- Lines 129-130: computes `report(spin_system,['1/max(diag(R))` using `report(spin_system,['1/max(diag(R)) = ' num2str(1/max_rate) ', tau_c = ' num2str(spin_system.rlx.tau_c{n})])`.
- Lines 130: computes `', tau_c` using `', tau_c = ' num2str(spin_system.rlx.tau_c{n})])`.
- Lines 146: computes `timer_nz` using `timer_nz=tic`.
- Lines 159: computes `rlx_shift` using `rlx_shift=sum(spin_system.chem.rp_rates)`.
- Lines 164: computes `error('nz_shift` using `error('nz_shift=''chem'' requires radical pair kinetics in inter.chem.')`.
- Lines 189: computes `omega_tau` using `omega_tau=cheap_norm(L0)*max(spin_system.rlx.tau_c{n})`.
- Lines 233: computes `Lp_left` using `Lp_left=operator(spin_system,{'L+'},{n},'left')`.
- Lines 234: computes `Lp_right` using `Lp_right=operator(spin_system,{'L+'},{n},'right')`.
- Lines 235: computes `Lm_left` using `Lm_left=operator(spin_system,{'L-'},{n},'left')`.
- Lines 236: computes `Lm_right` using `Lm_right=operator(spin_system,{'L-'},{n},'right')`.

### Local helper functions

- Line 797: `defaults()` — `function spin_system=defaults(spin_system)`.
  - Representative operation: `if ~isfield(spin_system.rlx,'equilibrium')`.
  - Representative operation: `report(spin_system,'relaxation destination not specified, assuming zero.')`.
- Line 815: `grumble()` — `function grumble(spin_system)`.
  - Representative operation: `if ~isfield(spin_system,'rlx')`.
  - Representative operation: `error('relaxation data (.rlx) is missing from the spin_system structure.')`.

## Parameters / inputs

- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- by those theories that support relaxation
- rate anisotropy. It has no effect on tho-
- se theories (e.g. Redfield) that do not.

## Outputs

- R -relaxation superoperator. If a Liouvillian is
- assembled manually, this dissipative superoperator
- must enter as 1i*R, for example
- L=H+1i*R+1i*K; do not use H+R+K.
- Note: a variety of relaxation theories are supported, see the relax-
- ation theory parameters section of the online manual.
- Note: Spinach context functions include relaxation and kinetics
- superoperators into the total Liovillian automatically.
- Note: dissipative bosonic modes declared in inter.modes receive
- thermalised GKSL dissipators from rlx_modes.m in Liouville
- space formalisms; the euler_angles parameter refers to the
- spin subsystem only and has no effect on the mode terms.

## Implementation structure

- Relaxation superoperator. Syntax:
- R=relaxation(spin_system,euler_angles)
- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- by those theories that support relaxation
- rate anisotropy. It has no effect on tho-
- se theories (e.g. Redfield) that do not.
- R -relaxation superoperator. If a Liouvillian is
- assembled manually, this dissipative superoperator
- must enter as 1i*R, for example
- L=H+1i*R+1i*K; do not use H+R+K.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `grumble()`, `mprealloc()`, `isfield()`, `any()`, `ismember()`, `report()`, `exist()`, `rlx_t1_t2()`, `theory()`, `hamiltonian()`, `assume()`, `num2str()`, `toc()`, `ischar()`, `point()`.
