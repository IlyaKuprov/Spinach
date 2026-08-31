# kernel/decouple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/decouple.m`
- Signature: `[L,rho]=decouple(spin_system,L,rho,spins)`
- Total lines: 225

## Purpose

Obliterates all interactions and populations in the subspace of states that involve the specified spins in any way. The specified spins would not contribute to the system dynamics until the Liouvillian is rebuilt from scratch. Syntax: [L,rho]=decouple(spin_system,L,rho,spins)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Return if the spin list is empty; implemented by `if isempty(spins), return; end`.
- Lines 51-52: Check consistency; implemented by `grumble(spin_system,L,rho,spins)`.
- Lines 54-55: Find the spins to be decoupled; implemented by `if isnumeric(spins)`.
- Lines 57-58: Find spins by numbers; implemented by `dec_mask=false(1,spin_system.comp.nspins); dec_mask(spins)=true`.
- Lines 62-63: Find spins by name; implemented by `dec_mask=ismember(spin_system.comp.isotopes,spins)`.
- Lines 67-68: Inform the user; implemented by `report(spin_system,[num2str(nnz(dec_mask)) ' spins to be frozen and depopulated.'])`.
- Lines 70-71: Get the spin space dimension; implemented by `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 73-74: Build the wipeout machinery; implemented by `switch spin_system.bas.formalism`.
- Lines 78-79: Get the list of states to be wiped; implemented by `zero_mask=(sum(spin_system.bas.basis(:,dec_mask),2)~=0)`.
- Lines 83-84: Hilbert space dimension and multiplicities; implemented by `dim=prod(spin_system.comp.mults)`.
- Lines 87-88: Identity component projector of the decoupled spins; implemented by `P=speye(dim^2)`.
- Lines 104-105: Process the Liouvillian; implemented by `if (nargout>0)&&(~isempty(L))`.
- Lines 107-108: Get dimension statistics; implemented by `spc_dim=size(L,1)/spn_dim`.
- Lines 110-111: Inform the user; implemented by `report(spin_system,['space sub-problem dimension: ' num2str(spc_dim)])`.
- Lines 114-115: Apply the wipeout machinery; implemented by `switch spin_system.bas.formalism`.
- Lines 119-120: Kron the list into the Fokker-Planck space; implemented by `fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask))`.
- Lines 122-124: Inform the user; implemented by `report(spin_system,['zeroing ' num2str(nnz(fp_zero_mask)) ' rows and columns in the Liouvillian.'])`.
- Lines 126-127: Apply the zero mask; implemented by `L(fp_zero_mask,:)=0; L(:,fp_zero_mask)=0`.

### Control flow inferred from the code

- Line 49: conditional branch on `isempty(spins), return; end`.
- Line 55: conditional branch on `isnumeric(spins)`.
- Line 74: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `{'zeeman-liouv','zeeman-hilb'}`, `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`, `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`.
- Line 89: `for` loop over `n=find(dec_mask)`.
- Line 91: `for` loop over `p=1:mults(n)`.
- Line 92: `for` loop over `q=1:mults(n)`.
- Line 105: conditional branch on `(nargout>0)&&(~isempty(L))`.
- Line 115: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`, `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`.
- Line 150: conditional branch on `(nargout>1)&&(~isempty(rho))`.
- Line 160: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'zeeman-hilb'`.

### Key state/data transformations

- Lines 58: computes `dec_mask` using `dec_mask=false(1,spin_system.comp.nspins); dec_mask(spins)=true`.
- Lines 71: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 84: computes `dim` using `dim=prod(spin_system.comp.mults)`.
- Lines 85: computes `mults` using `mults=spin_system.comp.mults`.
- Lines 88: computes `P` using `P=speye(dim^2)`.
- Lines 90: computes `idch` using `idch=sparse(dim^2,dim^2)`.
- Lines 93-95: computes `A` using `A=kron(kron(speye(prod(mults(1:(n-1)))), sparse(p,q,1,mults(n),mults(n))), speye(prod(mults((n+1):end))))`.
- Lines 108: computes `spc_dim` using `spc_dim=size(L,1)/spn_dim`.
- Lines 120: computes `fp_zero_mask` using `fp_zero_mask=logical(kron(ones(spc_dim,1),zero_mask))`.
- Lines 127: computes `L(fp_zero_mask,:)` using `L(fp_zero_mask,:)=0; L(:,fp_zero_mask)=0`.
- Lines 132: computes `fp_proj` using `fp_proj=kron(speye(spc_dim),P)`.
- Lines 135: computes `L` using `L=fp_proj*L*fp_proj`.
- Lines 172: computes `rho(fp_zero_mask,:)` using `rho(fp_zero_mask,:)=0`.
- Lines 180: computes `rho` using `rho=fp_proj*rho`.

### Local helper functions

- Line 194: `grumble()` — `function grumble(spin_system,L,rho,spins)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})`.
  - Representative operation: `error('analytical decoupling is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.')`.

## Parameters / inputs

- L -Liouvillian superoperator or, in zeeman-hilb,
- the Hamiltonian; this may be left empty
- rho -state vector or a horizontal stack thereof or,
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof; this may be left empty
- spins -spins to be wiped, specified either by name, e.g.
- {'13C','1H'}, or by a list of numbers, e.g. [1 2]

## Outputs

- rho -state vector(s) with all populations of the
- states involving the target spins set to zero
- L -Liouvillian superoperator with all terms in-
- volving the target spins set to zero
- Note: this function is an analytical equivalent of a perfect decoup-
- ling pulse sequence on the specified spins.
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported in
- the Liouville space formalisms. In the Zeeman formalisms the
- operation is an exact projection onto the subspace where the
- decoupled spins carry only their identity component, because
- spin involvement is not diagonal in the Zeeman basis. In
- zeeman-hilb, the Hamiltonian and the density matrices are
- stretched into Liouville space, projected there, and folded
- back; this replaces every decoupled-spin factor of the Hamil-
- tonian by its identity component average.

## Implementation structure

- Obliterates all interactions and populations in the subspace of states
- that involve the specified spins in any way. The specified spins would
- not contribute to the system dynamics until the Liouvillian is rebuilt
- from scratch. Syntax:
- [L,rho]=decouple(spin_system,L,rho,spins)
- L -Liouvillian superoperator or, in zeeman-hilb,
- the Hamiltonian; this may be left empty
- rho -state vector or a horizontal stack thereof or,
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof; this may be left empty
- spins -spins to be wiped, specified either by name, e.g.
- {'13C','1H'}, or by a list of numbers, e.g. [1 2]

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `false()`, `dec_mask()`, `ismember()`, `report()`, `num2str()`, `nnz()`, `speye()`, `mults()`, `conj()`, `logical()`, `clean_up()`, `rho()`, `iscell()`, `any()`.
