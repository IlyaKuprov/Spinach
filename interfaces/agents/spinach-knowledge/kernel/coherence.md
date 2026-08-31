# kernel/coherence.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/coherence.m`
- Signature: `rho=coherence(spin_system,rho,spec)`
- Total lines: 186

## Purpose

Coherence order selection function -keeps only the specified orders of coherence in the state vector. This is useful as an analytical re- placement for complicated phase cycles. Syntax: rho=coherence(spin_system,rho,spec)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(spin_system,rho,spec)`.
- Lines 46-47: Store dimension statistics; implemented by `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 54-55: Fold indirect dimensions; implemented by `rho=reshape(rho,[spn_dim spc_dim])`.
- Lines 57-58: Compute coherence order bookkeeping array; implemented by `switch spin_system.bas.formalism`.
- Lines 62-63: Projection quantum numbers of basis states; implemented by `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 67-68: Projection quantum numbers of ket and bra indices; implemented by `nspins=spin_system.comp.nspins`.
- Lines 73-74: Coherence orders of stretched density matrix elements; implemented by `M=M_ket-M_bra`.
- Lines 78-79: Projection quantum numbers of the Zeeman basis; implemented by `hdim=size(spin_system.bas.basis,1)`.
- Lines 83-84: Coherence orders of stretched density matrix elements; implemented by `M=repmat(M_lvl,[hdim 1])-kron(M_lvl,ones(hdim,1))`.
- Lines 88-89: Preallocate state mask array; implemented by `state_mask=false(spn_dim,numel(spec))`.
- Lines 91-92: Loop over specifications; implemented by `for n=1:numel(spec)`.
- Lines 94-95: Parse spin specification; implemented by `if ischar(spec{n}{1})`.
- Lines 97-98: Symbolic specification; implemented by `if strcmp(spec{n}{1},'all')`.
- Lines 110-111: Specification by number; implemented by `spins=spec{n}{1}`.
- Lines 115-116: Determine coherence order of each basis state; implemented by `coherence_orders_present=sum(M(:,spins),2)`.
- Lines 118-119: Wipe all coherence orders except those specified by the user; implemented by `state_mask(:,n)=ismember(coherence_orders_present,spec{n}{2})`.
- Lines 123-124: Intersect state masks; implemented by `state_mask=all(state_mask,2)`.
- Lines 126-127: Apply the state mask; implemented by `rho(~state_mask,:)=0`.

### Control flow inferred from the code

- Line 48: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 58: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`.
- Line 92: `for` loop over `n=1:numel(spec)`.
- Line 95: conditional branch on `ischar(spec{n}{1})`.
- Line 98: conditional branch on `strcmp(spec{n}{1},'all')`.
- Line 133: conditional branch on `norm(rho,1)<1e-10`.

### Key state/data transformations

- Lines 47: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 51: computes `spc_dim` using `spc_dim=numel(rho)/spn_dim`.
- Lines 52: computes `problem_dims` using `problem_dims=size(rho)`.
- Lines 55: computes `rho` using `rho=reshape(rho,[spn_dim spc_dim])`.
- Lines 63: computes `[~,M]` using `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 68: computes `nspins` using `nspins=spin_system.comp.nspins`.
- Lines 69: computes `spns` using `spns=(spin_system.comp.mults-1)/2`.
- Lines 70: computes `M_ket` using `M_ket=spns-spin_system.bas.basis(:,1:nspins)+1`.
- Lines 71: computes `M_bra` using `M_bra=spns-spin_system.bas.basis(:,(nspins+1):end)+1`.
- Lines 74: computes `M` using `M=M_ket-M_bra`.
- Lines 79: computes `hdim` using `hdim=size(spin_system.bas.basis,1)`.
- Lines 81: computes `M_lvl` using `M_lvl=spns-spin_system.bas.basis+1`.
- Lines 89: computes `state_mask` using `state_mask=false(spn_dim,numel(spec))`.
- Lines 99: computes `spins` using `spins=1:numel(spin_system.comp.isotopes)`.
- Lines 116: computes `coherence_orders_present` using `coherence_orders_present=sum(M(:,spins),2)`.
- Lines 119: computes `state_mask(:,n)` using `state_mask(:,n)=ismember(coherence_orders_present,spec{n}{2})`.
- Lines 127: computes `rho(~state_mask,:)` using `rho(~state_mask,:)=0`.

### Local helper functions

- Line 140: `grumble()` — `function grumble(spin_system,rho,spec)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})`.
  - Representative operation: `error('analytical coherence order selection is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.')`.

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- spec -a cell array containing the specification of
- which coherences to keep on which spins. For
- example
- {{'13C',[1 -1]},{'1H',-1}}
- keeps the states that have coherence order
- ((1 OR -1 on 13C) AND (-1 on 1H))
- instead of specific spins, it is possible to
- specify 'electrons', 'nuclei', and 'all'

## Outputs

- rho -the state vector with the undesired orders of
- spin correlations zeroed out
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported
- in the Liouville space formalisms. In zeeman-hilb, the densi-
- ty matrices are stretched into Liouville space, filtered the-
- re, and folded back.

## Implementation structure

- Coherence order selection function -keeps only the specified orders
- of coherence in the state vector. This is useful as an analytical re-
- placement for complicated phase cycles. Syntax:
- rho=coherence(spin_system,rho,spec)
- rho - a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- spec - a cell array containing the specification of
- which coherences to keep on which spins. For
- example
- {{'13C',[1 -1]},{'1H',-1}}
- keeps the states that have coherence order

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `lin2lm()`, `false()`, `ischar()`, `cellfun()`, `state_mask()`, `ismember()`, `all()`, `rho()`, `report()`, `vector()`, `iscell()`, `isvector()`, `any()`.
