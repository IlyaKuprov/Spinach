# kernel/hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/hamiltonian.m`
- Signature: `[I,Q]=hamiltonian(spin_system,operator_type)`
- Total lines: 1478

## Purpose

Hamiltonian operator or superoperator and its rotational decomposi- tion. Descriptor and operator generation are parallelised. Syntax: [I,Q]=hamiltonian(spin_system,operator_type)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `mode_quads()`, `spin_facts()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Set the default for the type; implemented by `if ~exist('operator_type','var'), operator_type='comm'; end`.
- Lines 49-50: Check consistency; implemented by `grumble(spin_system,operator_type)`.
- Lines 52-53: Decide if the Q part is required; implemented by `build_aniso=(nargout>1)`.
- Lines 55-56: Inform the user; implemented by `report(spin_system,'building Hamiltonian descriptor ')`.
- Lines 58-59: Preallocate spin number tables; implemented by `nL=zeros(spin_system.comp.nspins,8)`.
- Lines 62-63: Preallocate operator specifications; implemented by `opL(1:spin_system.comp.nspins,1:8)={'E'}`.
- Lines 66-67: Preallocate isotropic Hamiltonian coefficients; implemented by `isotropic=zeros(spin_system.comp.nspins,8)`.
- Lines 69-70: Preallocate spherical tensor coefficients; implemented by `ist_coeff{1}=zeros([spin_system.comp.nspins 8 3],'like',1i)`.
- Lines 73-74: Preallocate irreducible components; implemented by `irr_comp{1}=zeros([spin_system.comp.nspins 8 3],'like',1i)`.
- Lines 77-78: Process Zeeman interactions and NQI; implemented by `for n=1:spin_system.comp.nspins`.
- Lines 80-81: Write the isotropic Zeeman part; implemented by `switch spin_system.inter.zeeman.strength{n}`.
- Lines 85-86: Keep the carrier frequency; implemented by `zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 88-89: Update the Hamiltonian; implemented by `if abs(zeeman_iso)>spin_system.tols.liouv_zero`.
- Lines 91-92: Inform the user; implemented by `report(spin_system,['complete isotropic Zeeman interaction for spin ' num2str(n) ' '])`.
- Lines 95-96: Update the Hamiltonian descriptor; implemented by `nL(n,2)=n; opL(n,2)={'Lz'}; isotropic(n,2)=zeeman_iso`.
- Lines 102-104: Subtract the carrier frequency; implemented by `zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3- spin_system.inter.basefrqs(n)`.
- Lines 109-110: Inform the user; implemented by `report(spin_system,['offset isotropic Zeeman interaction for spin ' num2str(n) ' '])`.
- Lines 120-121: Inform the user; implemented by `report(spin_system,['isotropic Zeeman interaction ignored for spin ' num2str(n) '.'])`.

### Control flow inferred from the code

- Line 47: conditional branch on `~exist('operator_type','var'), operator_type='comm'; end`.
- Line 78: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 81: dispatches on `spin_system.inter.zeeman.strength{n}`; cases `'full'`, `'secular'`, `'ignore'`.
- Line 89: conditional branch on `abs(zeeman_iso)>spin_system.tols.liouv_zero`.
- Line 107: conditional branch on `abs(zeeman_iso)>spin_system.tols.liouv_zero`.
- Line 131: conditional branch on `build_aniso`.
- Line 137: conditional branch on `(norm(lam_zeeman,2)>spin_system.tols.liouv_zero)||`.
- Line 141: `for` loop over `k=1:3`.
- Line 147: dispatches on `spin_system.inter.zeeman.strength{n}`; cases `'full'`, `'secular'`, `'ignore'`.
- Line 153: conditional branch on `norm(lam_zeeman,2)>spin_system.tols.liouv_zero`.
- Line 157: conditional branch on `norm(phi_zeeman,2)>spin_system.tols.liouv_zero`.
- Line 192: conditional branch on `norm(spin_system.inter.coupling.matrix{n,n},2)>2*pi*spin_system.tols.inter_cutoff`.
- Line 198: conditional branch on `(norm(iso_quad,2)>1e-6*norm(phi_quad,2))||`.
- Line 204: `for` loop over `k=4:8`.

### Key state/data transformations

- Lines 53: computes `build_aniso` using `build_aniso=(nargout>1)`.
- Lines 59: computes `nL` using `nL=zeros(spin_system.comp.nspins,8)`.
- Lines 60: computes `nS` using `nS=zeros(spin_system.comp.nspins,8)`.
- Lines 63: computes `opL(1:spin_system.comp.nspins,1:8)` using `opL(1:spin_system.comp.nspins,1:8)={'E'}`.
- Lines 64: computes `opS(1:spin_system.comp.nspins,1:8)` using `opS(1:spin_system.comp.nspins,1:8)={'E'}`.
- Lines 67: computes `isotropic` using `isotropic=zeros(spin_system.comp.nspins,8)`.
- Lines 70: computes `ist_coeff{1}` using `ist_coeff{1}=zeros([spin_system.comp.nspins 8 3],'like',1i)`.
- Lines 71: computes `ist_coeff{2}` using `ist_coeff{2}=zeros([spin_system.comp.nspins 8 5],'like',1i)`.
- Lines 74: computes `irr_comp{1}` using `irr_comp{1}=zeros([spin_system.comp.nspins 8 3],'like',1i)`.
- Lines 75: computes `irr_comp{2}` using `irr_comp{2}=zeros([spin_system.comp.nspins 8 5],'like',1i)`.
- Lines 86: computes `zeeman_iso` using `zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 96: computes `nL(n,2)` using `nL(n,2)=n; opL(n,2)={'Lz'}; isotropic(n,2)=zeeman_iso`.
- Lines 134: computes `[~,lam_zeeman,phi_zeeman]` using `[~,lam_zeeman,phi_zeeman]=mat2sphten(spin_system.inter.zeeman.matrix{n})`.
- Lines 142: computes `irr_comp{1}(n,k,:)` using `irr_comp{1}(n,k,:)=lam_zeeman`.
- Lines 143: computes `irr_comp{2}(n,k,:)` using `irr_comp{2}(n,k,:)=phi_zeeman`.
- Lines 164: computes `nL(n,1)` using `nL(n,1)=n; opL(n,1)={'L+'}; ist_coeff{1}(n,1,1)=-0.5; ist_coeff{2}(n,1,2)=-0.5`.
- Lines 166: computes `nL(n,3)` using `nL(n,3)=n; opL(n,3)={'L-'}; ist_coeff{1}(n,3,3)=-0.5; ist_coeff{2}(n,3,4)=+0.5`.
- Lines 195: computes `[iso_quad,lam_quad,phi_quad]` using `[iso_quad,lam_quad,phi_quad]=mat2sphten(spin_system.inter.coupling.matrix{n,n})`.

### Local helper functions

- Line 663: `parfor_progr()` — `function parfor_progr()`. Build component operators in XYZ form
  - Representative operation: `terms_done=terms_done+1; last_message=toc-last_toc`.
  - Representative operation: `if (last_message>5)||(terms_done==nterms)`.
- Line 1384: `mode_quads()` — `function [mops,midx,mcoeff]=mode_quads(order,k1,k2)`. First derivative quadrature is sqrt(1/2)*(Cr+An)
  - Representative operation: `if order==1`.
  - Representative operation: `mops={{'C'},{'A'}}; midx={{k1},{k1}}`.
- Line 1409: `spin_facts()` — `function [sops,sidx,scoeff]=spin_facts(dblock,ns,ps)`. Effective field vectors expand into single-spin operators
  - Representative operation: `if isequal(size(dblock),[1 3])`.
  - Representative operation: `sops={{'L+'},{'L-'},{'Lz'}}; sidx={{ns},{ns},{ns}}`.
- Line 1444: `grumble()` — `function grumble(spin_system,operator_type)`. Catch operator type errors
  - Representative operation: `if ~isfield(spin_system,'bas')`.
  - Representative operation: `error('basis set information is missing, run basis() before calling this function.')`.

## Parameters / inputs

- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- in Hilbert space this parameter is ignored.

## Outputs

- I -rotationally invariant part of the Hamiltonian
- Q -irreducible components of the anisotropic part,
- use orientation.m to get the full Hamiltonian
- at each specific orientation
- Note: the code has a few rather eccentric blocks that bring the
- memory footprint to the absolute minimum and work around
- the sparse matrix addition efficiency problem.
- Note: bosonic mode terms declared in inter.modes are added to the
- invariant part I at their input orientation; the Q part and
- the orientation machinery refer to the spin subsystem only,
- for which rotations are well-defined.

## Implementation structure

- Hamiltonian operator or superoperator and its rotational decomposi-
- tion. Descriptor and operator generation are parallelised. Syntax:
- [I,Q]=hamiltonian(spin_system,operator_type)
- in Liouville space, operator_type can be set to
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- in Hilbert space this parameter is ignored.
- I -rotationally invariant part of the Hamiltonian
- Q -irreducible components of the anisotropic part,
- use orientation.m to get the full Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `report()`, `opL()`, `opS()`, `num2str()`, `isotropic()`, `mat2sphten()`, `lam_zeeman()`, `phi_zeeman()`, `phi_quad()`, `table()`, `clear()`, `cellfun()`, `pair_list()`, `lam_coupling()`.
