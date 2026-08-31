# kernel/conventions/transforms/ham2nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ham2nqi.m`
- Signature: `[omega,Q]=ham2nqi(H)`
- Total lines: 94

## Purpose

Converts a single-spin Hamiltonian back into the Zeeman and quadrupolar interaction parameters that had been used to generate it. Syntax: [omega,Q]=ham2nqi(H)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(H)`.
- Lines 38-39: Get multiplicity and basis operators; implemented by `mult=size(H,1); T=irr_sph_ten(mult); S=pauli(mult)`.
- Lines 41-42: Get Larmor frequencies (rad/s); implemented by `omega(1)=trace(S.x'*H)/norm(S.x,'fro')^2`.
- Lines 47-48: Spin > 1/2; implemented by `if mult>2`.
- Lines 50-51: Get quadratic expansion coefficients; implemented by `rank2(1)=trace(T{5}'*H)/norm(T{5},'fro')^2`.
- Lines 57-58: Translate quadratic part (rad/s); implemented by `Q=sphten2mat([],[],rank2); Q=real(Q)`.
- Lines 62-63: Spin 1/2; implemented by `Q=zeros(3,3)`.
- Lines 67-71: Run an explicit reconstruction; implemented by `HR=omega(1)*S.x+omega(2)*S.y+omega(3)*S.z+ Q(1,1)*S.x*S.x+Q(2,1)*S.y*S.x+Q(3,1)*S.z*S.x+ Q(1,2)*S.x*S.y+Q(2,2)*S.y*S.y+Q(3,2)*S.z*S.y+ Q(1,3)*S.x*S.z+Q(2,3)*S.y*S.z+Q(3,…`.
- Lines 73-74: Double-check reconstruction; implemented by `if norm(H-HR,2)>1e-6*norm(H,2)`.

### Control flow inferred from the code

- Line 48: conditional branch on `mult>2`.
- Line 74: conditional branch on `norm(H-HR,2)>1e-6*norm(H,2)`.

### Key state/data transformations

- Lines 39: computes `mult` using `mult=size(H,1); T=irr_sph_ten(mult); S=pauli(mult)`.
- Lines 42: computes `omega(1)` using `omega(1)=trace(S.x'*H)/norm(S.x,'fro')^2`.
- Lines 43: computes `omega(2)` using `omega(2)=trace(S.y'*H)/norm(S.y,'fro')^2`.
- Lines 44: computes `omega(3)` using `omega(3)=trace(S.z'*H)/norm(S.z,'fro')^2`.
- Lines 45: computes `omega` using `omega=real(omega)`.
- Lines 51: computes `rank2(1)` using `rank2(1)=trace(T{5}'*H)/norm(T{5},'fro')^2`.
- Lines 52: computes `rank2(2)` using `rank2(2)=trace(T{6}'*H)/norm(T{6},'fro')^2`.
- Lines 53: computes `rank2(3)` using `rank2(3)=trace(T{7}'*H)/norm(T{7},'fro')^2`.
- Lines 54: computes `rank2(4)` using `rank2(4)=trace(T{8}'*H)/norm(T{8},'fro')^2`.
- Lines 55: computes `rank2(5)` using `rank2(5)=trace(T{9}'*H)/norm(T{9},'fro')^2`.
- Lines 58: computes `Q` using `Q=sphten2mat([],[],rank2); Q=real(Q)`.
- Lines 68-71: computes `HR` using `HR=omega(1)*S.x+omega(2)*S.y+omega(3)*S.z+ Q(1,1)*S.x*S.x+Q(2,1)*S.y*S.x+Q(3,1)*S.z*S.x+ Q(1,2)*S.x*S.y+Q(2,2)*S.y*S.y+Q(3,2)*S.z*S.y+ Q(1,3)*S.x*S.z+Q(2,3)*S.y*S.z+Q(3,…`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(H)`. I made a mistake in a pulse program by typing "(360) 135135" instead of "(360) 135 135". I felt very guilty because I thought we'd have to
  - Representative operation: `if (~isnumeric(H))||(~ishermitian(H))|| (abs(trace(H))>eps()*norm(H,2))||(numel(H)<4)`.
  - Representative operation: `(abs(trace(H))>eps()*norm(H,2))||(numel(H)<4)`.

## Parameters / inputs

- H -single-spin Hamiltonian written in
- the Zeeman basis for a spin of any
- multiplicity

## Outputs

- omega -Larmor frequencies, rad/s
- Q -symmetric traceless quadrupolar
- coupling tensor, rad/s
- The outputs are returned such that:
- H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz +
- + [Sx Sy Sz]*Q*[Sx Sy Sz].';
- An error is produced if the Hamilonian contains
- any terms (for example, cubic) beyond those, or
- if it is not Hermitian and traceless.

## Implementation structure

- Converts a single-spin Hamiltonian back into the
- Zeeman and quadrupolar interaction parameters that
- had been used to generate it. Syntax:
- [omega,Q]=ham2nqi(H)
- H -single-spin Hamiltonian written in
- the Zeeman basis for a spin of any
- multiplicity
- omega -Larmor frequencies, rad/s
- Q -symmetric traceless quadrupolar
- coupling tensor, rad/s
- The outputs are returned such that:
- H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz +

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `irr_sph_ten()`, `pauli()`, `omega()`, `rank2()`, `sphten2mat()`, `ishermitian()`, `eps()`.
