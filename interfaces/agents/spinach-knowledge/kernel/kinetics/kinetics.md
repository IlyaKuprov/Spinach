# kernel/kinetics/kinetics.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/kinetics/kinetics.m`
- Signature: `K=kinetics(spin_system)`
- Total lines: 249

## Purpose

Chemical kinetics superoperator. Syntax: K=kinetics(spin_system)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Preallocate the answer; implemented by `K=mprealloc(spin_system,0)`.
- Lines 39-40: Find chemical reactions; implemented by `[sources,destins,rates]=find(spin_system.chem.rates)`.
- Lines 42-43: Loop over chemical reactions; implemented by `for n=1:numel(rates)`.
- Lines 45-46: Catch incorrect calls; implemented by `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Lines 50-51: Get the spins involved; implemented by `source_spins=spin_system.chem.parts{sources(n)}`.
- Lines 54-55: Get the states involved; implemented by `source_states=logical(sum(spin_system.bas.basis(:,source_spins),2))`.
- Lines 58-60: Make sure basis tables match; implemented by `if ~isequal(spin_system.bas.basis(source_states,source_spins), spin_system.bas.basis(destin_states,destin_spins))`.
- Lines 64-65: Move to integer indexing; implemented by `source_states=find(source_states)`.
- Lines 68-69: Update the kinetics superoperator; implemented by `K=K+rates(n)*sparse(source_states,destin_states,ones(size(source_states)),size(K,1),size(K,2))`.
- Lines 73-74: Find magnetization fluxes; implemented by `[source_spins,destin_spins,flux_rate]=find(spin_system.chem.flux_rate)`.
- Lines 76-77: Process magnetization fluxes; implemented by `if ~isempty(flux_rate)`.
- Lines 84-85: Index single-and multi-spin orders (sso and mso); implemented by `sso_state_mask=(sum(logical(spin_system.bas.basis),2)==1)`.
- Lines 88-89: Decide how to proceed; implemented by `switch spin_system.chem.flux_type`.
- Lines 98-99: Loop over the fluxes; implemented by `for n=1:numel(flux_rate)`.
- Lines 101-102: Find single-spin sources and destinations; implemented by `sso_source_state_mask=sso_state_mask&(spin_system.bas.basis(:,source_spins(n))~=0)`.
- Lines 105-106: Identify stationary states; implemented by `sso_static_state_mask=sso_source_state_mask&sso_destin_state_mask`.
- Lines 110-112: Make sure subspaces match; implemented by `if ~isequal(spin_system.bas.basis(sso_source_state_mask,source_spins(n)), spin_system.bas.basis(sso_destin_state_mask,destin_spins(n)))`.
- Lines 116-117: Move to integer indexing; implemented by `sso_source_state_index=find(sso_source_state_mask)`.

### Control flow inferred from the code

- Line 43: `for` loop over `n=1:numel(rates)`.
- Line 46: conditional branch on `~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 59: conditional branch on `~isequal(spin_system.bas.basis(source_states,source_spins),`.
- Line 77: conditional branch on `~isempty(flux_rate)`.
- Line 80: conditional branch on `~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 89: dispatches on `spin_system.chem.flux_type`; cases `'intramolecular'`, `'intermolecular'`.
- Line 99: `for` loop over `n=1:numel(flux_rate)`.
- Line 111: conditional branch on `~isequal(spin_system.bas.basis(sso_source_state_mask,source_spins(n)),`.
- Line 129: dispatches on `spin_system.chem.flux_type`; cases `'intramolecular'`, `'intermolecular'`.
- Line 143: conditional branch on `~isequal(spin_system.bas.basis(mso_source_state_mask,source_spins(n)),`.
- Line 185: conditional branch on `(~isempty(spin_system.chem.rp_theory))&&`.
- Line 189: conditional branch on `~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
- Line 208: dispatches on `spin_system.chem.rp_theory`; cases `'exponential'`, `'haberkorn'`, `'jones-hore'`.
- Line 231: conditional branch on `nnz(K)==0`.

### Key state/data transformations

- Lines 37: computes `K` using `K=mprealloc(spin_system,0)`.
- Lines 40: computes `[sources,destins,rates]` using `[sources,destins,rates]=find(spin_system.chem.rates)`.
- Lines 51: computes `source_spins` using `source_spins=spin_system.chem.parts{sources(n)}`.
- Lines 52: computes `destin_spins` using `destin_spins=spin_system.chem.parts{destins(n)}`.
- Lines 55: computes `source_states` using `source_states=logical(sum(spin_system.bas.basis(:,source_spins),2))`.
- Lines 56: computes `destin_states` using `destin_states=logical(sum(spin_system.bas.basis(:,destin_spins),2))`.
- Lines 74: computes `[source_spins,destin_spins,flux_rate]` using `[source_spins,destin_spins,flux_rate]=find(spin_system.chem.flux_rate)`.
- Lines 86: computes `mso_state_mask` using `mso_state_mask=(sum(logical(spin_system.bas.basis),2)>1 )`.
- Lines 106: computes `sso_static_state_mask` using `sso_static_state_mask=sso_source_state_mask&sso_destin_state_mask`.
- Lines 107: computes `sso_source_state_mask` using `sso_source_state_mask=xor(sso_source_state_mask,sso_static_state_mask)`.
- Lines 108: computes `sso_destin_state_mask` using `sso_destin_state_mask=xor(sso_destin_state_mask,sso_static_state_mask)`.
- Lines 117: computes `sso_source_state_index` using `sso_source_state_index=find(sso_source_state_mask)`.
- Lines 118: computes `sso_destin_state_index` using `sso_destin_state_index=find(sso_destin_state_mask)`.
- Lines 138: computes `mso_static_state_mask` using `mso_static_state_mask=mso_source_state_mask&mso_destin_state_mask`.
- Lines 139: computes `mso_source_state_mask` using `mso_source_state_mask=xor(mso_source_state_mask,mso_static_state_mask)`.
- Lines 140: computes `mso_destin_state_mask` using `mso_destin_state_mask=xor(mso_destin_state_mask,mso_static_state_mask)`.
- Lines 149: computes `mso_source_state_index` using `mso_source_state_index=find(mso_source_state_mask)`.
- Lines 150: computes `mso_destin_state_index` using `mso_destin_state_index=find(mso_destin_state_mask)`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual. All adjustable pa-
- rameters are described in the chemical
- kinetics parameters section.

## Outputs

- K -kinetics superoperator. If a Liouvillian is
- assembled manually, this dissipative super-
- operator must enter as 1i*K, for example
- L=H+1i*R+1i*K
- Note: a large variety of chemical reaction models is supported,
- see the chemical kinetics parameters section of the onli-
- ne manual.
- Note: Spinach context functions include relaxation and kinetics
- superoperators into the total Liovillian automatically.

## Implementation structure

- Chemical kinetics superoperator. Syntax:
- K=kinetics(spin_system)
- spin_system - Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual. All adjustable pa-
- rameters are described in the chemical
- kinetics parameters section.
- K - kinetics superoperator. If a Liouvillian is
- assembled manually, this dissipative super-
- operator must enter as 1i*K, for example
- L=H+1i*R+1i*K

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mprealloc()`, `strcmp()`, `sources()`, `destins()`, `logical()`, `isequal()`, `rates()`, `report()`, `source_spins()`, `destin_spins()`, `xor()`, `flux_rate()`, `ismember()`, `unit_oper()`, `operator()`, `num2cell()`.
