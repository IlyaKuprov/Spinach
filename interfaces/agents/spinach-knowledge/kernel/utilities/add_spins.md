# kernel/utilities/add_spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/add_spins.m`
- Signature: `[mult,proj]=add_spins(spin_a,spin_b)`
- Total lines: 106

## Purpose

Reduction of direct products of two su(2) irreps. Syntax: [mult,proj]=add_spins(spin_a,spin_b)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_a,spin_b)`.
- Lines 32-33: Get the individual irreps; implemented by `spin_a=pauli(2*spin_a+1)`.
- Lines 36-37: Build direct product rep generators; implemented by `Sx=kron(spin_a.u,spin_b.x)+kron(spin_a.x,spin_b.u)`.
- Lines 41-42: Diagonalise Casimir operator; implemented by `[V,D]=eig(full(Sx^2+Sy^2+Sz^2),'vector')`.
- Lines 44-45: Index Casimir operator eigenvalues; implemented by `[S_sq,~,idx_b]=unique(uint32(D))`.
- Lines 47-48: Make projectors; implemented by `proj=cell(1,numel(S_sq))`.
- Lines 53-54: Canonicalise projectors; implemented by `for n=1:numel(S_sq)`.
- Lines 56-57: Set the expectations; implemented by `mult(n)=size(proj{n},2)`.
- Lines 60-61: Set Sz eigenvalue order; implemented by `Sz_irr=proj{n}'*Sz*proj{n}`.
- Lines 65-66: Require diagonal Sz; implemented by `V=V(:,idx); proj{n}=proj{n}*V`.
- Lines 71-72: Require real and positive Sx; implemented by `while any(proj{n}'*Sx*proj{n}+sqrt(eps)<0,'all')`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:numel(S_sq)`.
- Line 54: `for` loop over `n=1:numel(S_sq)`.
- Line 67: conditional branch on `norm(proj{n}'*Sz*proj{n}-Q.z,1)>sqrt(eps)`.
- Line 72: `while` loop over `any(proj{n}'*Sx*proj{n}+sqrt(eps)<0,'all')`.
- Line 74: `for` loop over `k=1:size(Sx_irr,1)`.
- Line 75: conditional branch on `any(Sx_irr(:,k)+sqrt(eps)<0)`.
- Line 80: conditional branch on `(norm(proj{n}'*Sx*proj{n}-Q.x,1)>sqrt(eps))||`.

### Key state/data transformations

- Lines 33: computes `spin_a` using `spin_a=pauli(2*spin_a+1)`.
- Lines 34: computes `spin_b` using `spin_b=pauli(2*spin_b+1)`.
- Lines 37: computes `Sx` using `Sx=kron(spin_a.u,spin_b.x)+kron(spin_a.x,spin_b.u)`.
- Lines 38: computes `Sy` using `Sy=kron(spin_a.u,spin_b.y)+kron(spin_a.y,spin_b.u)`.
- Lines 39: computes `Sz` using `Sz=kron(spin_a.u,spin_b.z)+kron(spin_a.z,spin_b.u)`.
- Lines 42: computes `[V,D]` using `[V,D]=eig(full(Sx^2+Sy^2+Sz^2),'vector')`.
- Lines 45: computes `[S_sq,~,idx_b]` using `[S_sq,~,idx_b]=unique(uint32(D))`.
- Lines 48: computes `proj` using `proj=cell(1,numel(S_sq))`.
- Lines 57: computes `mult(n)` using `mult(n)=size(proj{n},2)`.
- Lines 58: computes `Q` using `Q=pauli(mult(n))`.
- Lines 61: computes `Sz_irr` using `Sz_irr=proj{n}'*Sz*proj{n}`.
- Lines 63: computes `[~,idx]` using `[~,idx]=sort(D,'descend')`.
- Lines 66: computes `V` using `V=V(:,idx); proj{n}=proj{n}*V`.
- Lines 73: computes `Sx_irr` using `Sx_irr=proj{n}'*Sx*proj{n}`.
- Lines 76: computes `proj{n}(:,k)` using `proj{n}(:,k)=-proj{n}(:,k); break`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(spin_a,spin_b)`.
  - Representative operation: `if (~isnumeric(spin_a))||(~isreal(spin_a))|| (~isscalar(spin_a))||(mod(2*spin_a+1,1)~=0)|| (spin_a<1/2)`.
  - Representative operation: `(~isscalar(spin_a))||(mod(2*spin_a+1,1)~=0)|| (spin_a<1/2)`.

## Parameters / inputs

- spin_a -quantum number of the first spin,
- an integer or a half-integer
- spin_b -quantum number of the second spin,
- an integer or a half-integer

## Outputs

- mult -multiplicities corresponding to
- the values of the total spin that
- are present
- proj -projectors that reduce the direct
- product representation, a cell ar-
- ray of matrices

## Implementation structure

- Reduction of direct products of two su(2) irreps. Syntax:
- [mult,proj]=add_spins(spin_a,spin_b)
- spin_a -quantum number of the first spin,
- an integer or a half-integer
- spin_b -quantum number of the second spin,
- mult -multiplicities corresponding to
- the values of the total spin that
- are present
- proj -projectors that reduce the direct
- product representation, a cell ar-
- ray of matrices
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `uint32()`, `mult()`, `any()`, `Sx_irr()`, `isscalar()`.
