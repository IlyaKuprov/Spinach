# experiments/pseudocon/geffect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/geffect.m`
- Signature: `g=geffect(spin_system,states)`
- Total lines: 125

## Purpose

Effective g-tensor for the user-specified Kramers doublet, computed as described in

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(spin_system,states)`.
- Lines 31-32: Get the g-tensor for each spin; implemented by `for n=spin_system.comp.nspins:-1:1`.
- Lines 36-37: Get Sx, Sy, Sz operators for each spin; implemented by `for n=spin_system.comp.nspins:-1:1`.
- Lines 45-46: Get magnetic moment operators; implemented by `Mx=sparse(0); My=sparse(0); Mz=sparse(0)`.
- Lines 53-54: Get the complete Hamiltonian; implemented by `[H,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 57-58: Find the eigenstates; implemented by `[V,D]=eig(full(H)); D=diag(D)`.
- Lines 60-61: Sort the eigenstates; implemented by `[~,index]=sort(D,'ascend'); V=V(:,index)`.
- Lines 63-64: Pick out the states; implemented by `V=V(:,states)`.
- Lines 77-78: Matrix square root; implemented by `g=real(sqrtm(G))`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=spin_system.comp.nspins:-1:1`.
- Line 37: `for` loop over `n=spin_system.comp.nspins:-1:1`.
- Line 47: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 70: `for` loop over `n=1:3`.
- Line 71: `for` loop over `k=1:3`.

### Key state/data transformations

- Lines 33: computes `g{n}` using `g{n}=gtensorof(spin_system,n)`.
- Lines 38-39: computes `Sx{n}` using `Sx{n}=(operator(spin_system,{'L+'},{n})+ operator(spin_system,{'L-'},{n}))/2`.
- Lines 40-41: computes `Sy{n}` using `Sy{n}=(operator(spin_system,{'L+'},{n})- operator(spin_system,{'L-'},{n}))/2i`.
- Lines 42: computes `Sz{n}` using `Sz{n}=operator(spin_system,{'Lz'},{n})`.
- Lines 46: computes `Mx` using `Mx=sparse(0); My=sparse(0); Mz=sparse(0)`.
- Lines 49: computes `My` using `My=My-g{n}(2,1)*Sx{n}-g{n}(2,2)*Sy{n}-g{n}(2,3)*Sz{n}`.
- Lines 50: computes `Mz` using `Mz=Mz-g{n}(3,1)*Sx{n}-g{n}(3,2)*Sy{n}-g{n}(3,3)*Sz{n}`.
- Lines 54: computes `[H,Q]` using `[H,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 55: computes `H` using `H=H+orientation(Q,[0 0 0]); H=(H+H')/2`.
- Lines 58: computes `[V,D]` using `[V,D]=eig(full(H)); D=diag(D)`.
- Lines 61: computes `[~,index]` using `[~,index]=sort(D,'ascend'); V=V(:,index)`.
- Lines 64: computes `V` using `V=V(:,states)`.
- Lines 67: computes `phi{1}` using `phi{1}=V'*Mx*V; phi{2}=V'*My*V; phi{3}=V'*Mz*V`.
- Lines 72-73: computes `G(n,k)` using `G(n,k)=(phi{n}(1,1)-phi{n}(2,2))*(phi{k}(1,1)-phi{k}(2,2))+ 2*phi{n}(1,2)*phi{k}(2,1)+2*phi{k}(1,2)*phi{n}(2,1)`.
- Lines 78: computes `g` using `g=real(sqrtm(G))`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(spin_system,states)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
  - Representative operation: `error('this function is only available in Hilbert space.')`.

## Syntax

```matlab
g=geffect(spin_system,states)
```

## Parameters / inputs

- states -the numbers of the states to use
- (numbered sequentially from the
- lowest to the highest energy)

## Outputs

- g -3x3 g-tensor matrix in Bohr mag-
- neton units

## Implementation structure

- Effective g-tensor for the user-specified Kramers
- doublet, computed as described in
- g=geffect(spin_system,states)
- states -the numbers of the states to use
- (numbered sequentially from the
- lowest to the highest energy)
- g -3x3 g-tensor matrix in Bohr mag-
- neton units
- Check consistency
- Get the g-tensor for each spin
- Get Sx, Sy, Sz operators for each spin
- Get magnetic moment operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `orientation()`, `sqrtm()`, `strcmp()`, `any()`, `states()`.
