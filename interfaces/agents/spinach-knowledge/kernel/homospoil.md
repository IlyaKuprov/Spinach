# kernel/homospoil.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/homospoil.m`
- Signature: `rho=homospoil(spin_system,rho,zqc_flag)`
- Total lines: 118

## Purpose

Emulates a strong homospoil pulse -only zero-frequency states with respect to the carrier frequencies (chemical shifts are not conside- red) survive the process. Syntax: rho=homospoil(spin_system,rho,zqc_flag)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(rho,zqc_flag)`.
- Lines 46-47: In Hilbert space, only keep the diagonal; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Lines 52-53: In Zeeman basis of Liouville space, zero off-diagonals in every stacked column; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Lines 60-61: Store dimension statistics; implemented by `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 65-66: Fold indirect dimensions; implemented by `rho=reshape(rho,[spn_dim spc_dim])`.
- Lines 68-69: Pull the projection information from the basis; implemented by `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 71-72: Filter the state vector; implemented by `switch zqc_flag`.
- Lines 76-77: Find the states that have zero carrier frequency and kill everything else; implemented by `rho(abs(sum(repmat(spin_system.inter.basefrqs,size(spin_system.bas.basis,1),1).*M,2))>1e-6,:)=0`.
- Lines 81-82: Find the longitudinal states and kill everything else; implemented by `rho(sum(abs(M),2)>0,:)=0`.
- Lines 86-87: Complain and bomb out; implemented by `error('unknown ZQC flag.')`.
- Lines 91-92: Unfold indirect dimensions; implemented by `rho=reshape(rho,problem_dims)`.
- Lines 94-95: Report overly destructive calls; implemented by `if norm(rho,1)<1e-10`.

### Control flow inferred from the code

- Line 47: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 53: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Line 72: dispatches on `zqc_flag`; cases `'keep'`, `'destroy'`.
- Line 95: conditional branch on `norm(rho,1)<1e-10`.

### Key state/data transformations

- Lines 48: computes `dim` using `dim=size(rho,1); rhod=diag(rho)`.
- Lines 49: computes `rho` using `rho=spdiags(rhod,0,dim,dim); return`.
- Lines 55: computes `offdiag_mask` using `offdiag_mask=true(dim^2,1)`.
- Lines 56: computes `offdiag_mask(1:(dim+1):dim^2)` using `offdiag_mask(1:(dim+1):dim^2)=false`.
- Lines 57: computes `rho(offdiag_mask,:)` using `rho(offdiag_mask,:)=0; return`.
- Lines 61: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 62: computes `spc_dim` using `spc_dim=numel(rho)/spn_dim`.
- Lines 63: computes `problem_dims` using `problem_dims=size(rho)`.
- Lines 69: computes `[~,M]` using `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 82: computes `rho(sum(abs(M),2)>0,:)` using `rho(sum(abs(M),2)>0,:)=0`.

### Local helper functions

- Line 102: `grumble()` — `function grumble(rho,zqc_flag)`.
  - Representative operation: `if ~isnumeric(rho)`.
  - Representative operation: `error('the state vector(s) must be numeric.')`.

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof
- zqc_flag -a flag controlling the fate of zero-quantum
- coherences. If set to 'keep', causes ZQCs to
- survive the process, approximating experimen-
- tal behaviour. If set to 'destroy', wipes the
- zero-quantum coherences -only the longitudi-
- nal states survive the process.
- The flag is ignored in zeeman-hilb and zeeman-
- liouv formalisms, where the effect is always
- to destroy everything except the diagonal of
- the density matrix.

## Outputs

- rho -the state vector(s) with only the longitudi-
- nal or only the zero-quantum states kept
- Note: this function is only available for sphten-liouv formalism; it
- supports Fokker-Planck direct products.
- Note: this is a purely mathematical filter that only mimics -in an
- idealised way -the effect of a real homospoil pulse. Essenti-
- ally, it searches the density matrix for any transverse state
- populations and zeroes them out. If the flag is set, zero-qua-
- ntum coherences are also erased.

## Implementation structure

- Emulates a strong homospoil pulse -only zero-frequency states with
- respect to the carrier frequencies (chemical shifts are not conside-
- red) survive the process. Syntax:
- rho=homospoil(spin_system,rho,zqc_flag)
- rho -a state vector or a horizontal stack thereof
- zqc_flag -a flag controlling the fate of zero-quantum
- coherences. If set to 'keep', causes ZQCs to
- survive the process, approximating experimen-
- tal behaviour. If set to 'destroy', wipes the
- zero-quantum coherences -only the longitudi-
- nal states survive the process.
- The flag is ignored in zeeman-hilb and zeeman-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `spdiags()`, `true()`, `offdiag_mask()`, `rho()`, `lin2lm()`, `report()`, `vector()`, `ischar()`, `ismember()`.
