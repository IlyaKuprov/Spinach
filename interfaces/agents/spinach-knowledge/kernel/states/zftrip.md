# kernel/states/zftrip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/zftrip.m`
- Signature: `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)`
- Total lines: 122

## Purpose

Projection of the zero-field triplet state with user-specified populations of Cartesian ZFS eigenstates onto the higher-field ZFS + Zeeman eigenstates. This is commonly seen in triplet DNP with photo-generated two-electron triplets. Syntax: rho=zftrip(spin_system,ZFS,pops,Z,B,idx)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(spin_system,ZFS,pops,Z,B,idx)`.
- Lines 44-45: Get electron spin-1 operators; implemented by `A=pauli(3); Ex=A.x; Ey=A.y; Ez=A.z`.
- Lines 47-50: Build ZFS Hamiltonian in the laboratory frame; implemented by `H_ZFS=2*pi*(ZFS(1,1)*Ex*Ex+ZFS(1,2)*Ex*Ey+ZFS(1,3)*Ex*Ez+ ZFS(2,1)*Ey*Ex+ZFS(2,2)*Ey*Ey+ZFS(2,3)*Ey*Ez+ ZFS(3,1)*Ez*Ex+ZFS(3,2)*Ez*Ey+ZFS(3,3)*Ez*Ez)`.
- Lines 52-53: Diagonalise ZFS Hamiltonian; implemented by `[V,D]=eig(full(H_ZFS),'vector')`.
- Lines 55-56: Match literature convention (|zz|>|xx|>|yy|); implemented by `[~,ord]=sort(abs(D),'ascend'); V=V(:,ord)`.
- Lines 59-60: Build zero field density matrix from populations; implemented by `DM=Vx*pops(1)*Vx'+Vy*pops(2)*Vy'+Vz*pops(3)*Vz'`.
- Lines 62-63: Build Zeeman Hamiltonian in the laboratory frame; implemented by `Bx=0; By=0; Bz=B`.
- Lines 68-69: Diagonalise high field Hamiltonian; implemented by `[V,~]=eig(full(H_Z+H_ZFS),'vector')`.
- Lines 71-72: Drop high-field coherences; implemented by `DM=V*diag(diag(V'*DM*V))*V'`.
- Lines 74-75: Construct a Spinach state; implemented by `T=irr_sph_ten(3); rho=complex(0)`.

### Control flow inferred from the code

- Line 76: `for` loop over `n=1:numel(T), T{n}=T{n}/norm(T{n},'fro')^2; end`.

### Key state/data transformations

- Lines 45: computes `A` using `A=pauli(3); Ex=A.x; Ey=A.y; Ez=A.z`.
- Lines 48-50: computes `H_ZFS` using `H_ZFS=2*pi*(ZFS(1,1)*Ex*Ex+ZFS(1,2)*Ex*Ey+ZFS(1,3)*Ex*Ez+ ZFS(2,1)*Ey*Ex+ZFS(2,2)*Ey*Ey+ZFS(2,3)*Ey*Ez+ ZFS(3,1)*Ez*Ex+ZFS(3,2)*Ez*Ey+ZFS(3,3)*Ez*Ez)`.
- Lines 53: computes `[V,D]` using `[V,D]=eig(full(H_ZFS),'vector')`.
- Lines 56: computes `[~,ord]` using `[~,ord]=sort(abs(D),'ascend'); V=V(:,ord)`.
- Lines 57: computes `Vz` using `Vz=V(:,3); Vx=V(:,2); Vy=V(:,1)`.
- Lines 60: computes `DM` using `DM=Vx*pops(1)*Vx'+Vy*pops(2)*Vy'+Vz*pops(3)*Vz'`.
- Lines 63: computes `Bx` using `Bx=0; By=0; Bz=B`.
- Lines 64-66: computes `H_Z` using `H_Z=2*pi*(Z(1,1)*Ex*Bx+Z(1,2)*Ex*By+Z(1,3)*Ex*Bz+ Z(2,1)*Ey*Bx+Z(2,2)*Ey*By+Z(2,3)*Ey*Bz+ Z(3,1)*Ez*Bx+Z(3,2)*Ez*By+Z(3,3)*Ez*Bz)`.
- Lines 69: computes `[V,~]` using `[V,~]=eig(full(H_Z+H_ZFS),'vector')`.
- Lines 75: computes `T` using `T=irr_sph_ten(3); rho=complex(0)`.
- Lines 77: computes `rho` using `rho=rho+trace(T{1}'*DM)*state(spin_system,'T0,0' ,idx)`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(spin_system,ZFS,pops,Z,B,idx)`.
  - Representative operation: `if (~isnumeric(ZFS))||(~isreal(ZFS))||(size(ZFS,1)~=3)|| (size(ZFS,2)~=3)||(norm(ZFS-ZFS','fro')>1e-6*norm(ZFS,'fro'))`.
  - Representative operation: `(size(ZFS,2)~=3)||(norm(ZFS-ZFS','fro')>1e-6*norm(ZFS,'fro'))`.

## Parameters / inputs

- ZFS -3x3 ZFS tensor (Hz) in the laboratory frame
- of reference; use zfs2mat() to get it from
- D, E, and molecular Euler angles
- pops -a three-element vector with populations of
- X, Y, and Z eigenstates of the ZFS tensor
- at zero magnetic field, order: [pX pY pZ]
- Z -3x3 Zeeman interaction tensor (Hz/Tesla) in
- the laboratory frame of reference; use func-
- tions like axrh2mat() to get it from eigen-
- values and molecular Euler angles
- B -magnetic field directed along the Z axis of
- the laboratory frame of reference, Tesla
- idx -index of the electron triplet (use 'E3') in
- the sys.isotopes list

## Outputs

- rho -spin density matrix (Hilbert space) or state
- vector (Liouville space)

## Implementation structure

- Projection of the zero-field triplet state with user-specified
- populations of Cartesian ZFS eigenstates onto the higher-field
- ZFS + Zeeman eigenstates. This is commonly seen in triplet DNP
- with photo-generated two-electron triplets. Syntax:
- rho=zftrip(spin_system,ZFS,pops,Z,B,idx)
- ZFS -3x3 ZFS tensor (Hz) in the laboratory frame
- of reference; use zfs2mat() to get it from
- D, E, and molecular Euler angles
- pops -a three-element vector with populations of
- X, Y, and Z eigenstates of the ZFS tensor
- at zero magnetic field, order: [pX pY pZ]
- Z -3x3 Zeeman interaction tensor (Hz/Tesla) in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `ZFS()`, `pops()`, `irr_sph_ten()`, `complex()`, `state()`, `any()`, `isscalar()`, `strcmp()`, `int2str()`.
