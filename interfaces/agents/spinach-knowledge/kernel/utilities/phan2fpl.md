# kernel/utilities/phan2fpl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/phan2fpl.m`
- Signature: `rho=phan2fpl(phan,rho)`
- Total lines: 48

## Purpose

Projects a spatial intensity distribution into the Fokker-Planck space, using it as the image painted by the the spin state supp- lied. Syntax: rho=phan2fpl(phan,rho)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(phan,rho)`.
- Lines 27-28: Stretch the phantom and kron it with the spin state; implemented by `rho=kron(phan(:),rho)`.

### Key state/data transformations

- Lines 28: computes `rho` using `rho=kron(phan(:),rho)`.

### Local helper functions

- Line 33: `grumble()` — `function grumble(phan,rho)`. Q: "How many members of a certain demographic group does it take to perform a specified task?"
  - Representative operation: `if (~isnumeric(rho))||(size(rho,2)~=1)`.
  - Representative operation: `error('rho must be a column vector.')`.

## Parameters / inputs

- phan -phantom (the spatial distribution of the
- amplitude of the specified spin state)
- rho -Liouville space state vector

## Outputs

- rho -Fokker-Planck state vector

## Implementation structure

- Projects a spatial intensity distribution into the Fokker-Planck
- space, using it as the image painted by the the spin state supp-
- lied. Syntax:
- rho=phan2fpl(phan,rho)
- phan -phantom (the spatial distribution of the
- amplitude of the specified spin state)
- rho -Liouville space state vector
- rho -Fokker-Planck state vector
- Check consistency
- Stretch the phantom and kron it with the spin state
- Consistency enforcement
- Q: "How many members of a certain demographic group does

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phan()`, `ismember()`, `ndims()`.
