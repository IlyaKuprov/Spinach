# kernel/utilities/sphten2zeeman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sphten2zeeman.m`
- Signature: `P=sphten2zeeman(spin_system)`
- Total lines: 77

## Purpose

Returns a matrix that converts state vectors written in the spherical tensor basis set used by Spinach into state vectors written in the Zeeman basis set in Liouville space. Syntax: P=sphten2zeeman(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system)`.
- Lines 32-34: Preallocate the answer; implemented by `P=spalloc(prod(spin_system.comp.mults.^2), size(spin_system.bas.basis,1),0)`.
- Lines 36-37: Destination basis is not normalised; implemented by `destin_norm=sqrt(prod(spin_system.comp.mults))`.
- Lines 39-40: Loop over the basis set; implemented by `parfor n=1:size(spin_system.bas.basis,1)`.
- Lines 42-43: Get the state going; implemented by `rho=sparse(1)`.
- Lines 45-46: Loop over the elements; implemented by `for k=1:size(spin_system.bas.basis,2)`.
- Lines 48-49: Get the spherical tensors for the current spin; implemented by `ists=irr_sph_ten(spin_system.comp.mults(k))`.
- Lines 51-52: Multiply into the state; implemented by `rho=kron(rho,ists{spin_system.bas.basis(n,k)+1})`.
- Lines 56-57: Source basis is not normalised; implemented by `source_norm=norm(rho(:),2)`.
- Lines 59-60: Write a column into the projector; implemented by `P(:,n)=destin_norm*rho(:)/source_norm`.

### Control flow inferred from the code

- Line 40: `parfor` loop over `n=1:size(spin_system.bas.basis,1)`.
- Line 46: `for` loop over `k=1:size(spin_system.bas.basis,2)`.

### Key state/data transformations

- Lines 33-34: computes `P` using `P=spalloc(prod(spin_system.comp.mults.^2), size(spin_system.bas.basis,1),0)`.
- Lines 37: computes `destin_norm` using `destin_norm=sqrt(prod(spin_system.comp.mults))`.
- Lines 43: computes `rho` using `rho=sparse(1)`.
- Lines 49: computes `ists` using `ists=irr_sph_ten(spin_system.comp.mults(k))`.
- Lines 57: computes `source_norm` using `source_norm=norm(rho(:),2)`.
- Lines 60: computes `P(:,n)` using `P(:,n)=destin_norm*rho(:)/source_norm`.

### Local helper functions

- Line 67: `grumble()` — `function grumble(spin_system)`. Nurture your minds with great thoughts. To believe in the heroic makes heroes.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- spin_system -main Spinach data structure using
- sphten-liouv formalism and inclu-
- ding basis set information

## Outputs

- P -projector matrix that is to be used in the fol-
- lowing way:
- rho_zeeman=P*rho_sphten
- Note: the projector need not be square and may be huge.

## Implementation structure

- Returns a matrix that converts state vectors written in the
- spherical tensor basis set used by Spinach into state vectors
- written in the Zeeman basis set in Liouville space. Syntax:
- P=sphten2zeeman(spin_system)
- spin_system -main Spinach data structure using
- sphten-liouv formalism and inclu-
- ding basis set information
- P -projector matrix that is to be used in the fol-
- lowing way:
- rho_zeeman=P*rho_sphten
- Note: the projector need not be square and may be huge.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `irr_sph_ten()`, `rho()`, `strcmp()`.
