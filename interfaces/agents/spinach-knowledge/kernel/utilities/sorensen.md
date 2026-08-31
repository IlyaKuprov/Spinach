# kernel/utilities/sorensen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sorensen.m`
- Signature: `b=sorensen(rho_init,rho_targ)`
- Total lines: 62

## Purpose

Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators. Equation 186 from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax: b=sorensen(rho_init,rho_targ)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(rho_init,rho_targ)`.
- Lines 31-32: Diagonalise both matrices; implemented by `[~,sigma_init]=eig(full(rho_init))`.
- Lines 35-36: Sort eigenvalues; implemented by `sigma_init=sort(diag(sigma_init))`.
- Lines 39-40: Compute Sorensen bound; implemented by `b=(sigma_init'*sigma_targ)/trace(rho_targ^2)`.

### Key state/data transformations

- Lines 32: computes `[~,sigma_init]` using `[~,sigma_init]=eig(full(rho_init))`.
- Lines 33: computes `[~,sigma_targ]` using `[~,sigma_targ]=eig(full(rho_targ))`.
- Lines 36: computes `sigma_init` using `sigma_init=sort(diag(sigma_init))`.
- Lines 37: computes `sigma_targ` using `sigma_targ=sort(diag(sigma_targ))`.
- Lines 40: computes `b` using `b=(sigma_init'*sigma_targ)/trace(rho_targ^2)`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(rho_init,rho_targ)`. People think they do not understand mathematics, but that rather depends on how you explain. If you ask an alcoholic which is big-
  - Representative operation: `if (~isnumeric(rho_init))||(~isnumeric(rho_targ))|| (~ishermitian(rho_init))||(~ishermitian(rho_targ))|| (numel(rho_init)~=numel(rho_targ))`.
  - Representative operation: `(~ishermitian(rho_init))||(~ishermitian(rho_targ))|| (numel(rho_init)~=numel(rho_targ))`.

## Parameters / inputs

- rho_init -initia ldensity matrix, Hilbert space
- rho_targ -target density matrix, Hilbert space
- Output:
- b -Sorensen bound
- Note: this is an exact unitary bound; the amount reachable
- with realistically available instrumental controls
- may be smaller, see the detailed analysis here:

## Implementation structure

- Sorensen bound for the maximum transfer efficiency between
- two states under arbitrary control operators. Equation 186
- from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax:
- b=sorensen(rho_init,rho_targ)
- rho_init -initia ldensity matrix, Hilbert space
- rho_targ -target density matrix, Hilbert space
- Output:
- b -Sorensen bound
- Note: this is an exact unitary bound; the amount reachable
- with realistically available instrumental controls
- may be smaller, see the detailed analysis here:
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ishermitian()`.
