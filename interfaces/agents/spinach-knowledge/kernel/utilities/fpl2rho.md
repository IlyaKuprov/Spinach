# kernel/utilities/fpl2rho.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/fpl2rho.m`
- Signature: `rho=fpl2rho(rho,dims)`
- Total lines: 55

## Purpose

Returns the average of the spin state vector across the spatial dimensions of the sample. Syntax: rho=fpl2rho(rho,dims)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(rho,dims)`.
- Lines 27-28: Find out the stack size; implemented by `stack_size=size(rho,2)`.
- Lines 30-31: Expose the spin dimension (no ND sparse support yet); implemented by `rho=reshape(full(rho),[size(rho,1)/prod(dims) prod(dims) stack_size])`.
- Lines 33-34: Average over the spatial coordinates; implemented by `rho=squeeze(sum(rho,2)/prod(dims))`.

### Key state/data transformations

- Lines 28: computes `stack_size` using `stack_size=size(rho,2)`.
- Lines 31: computes `rho` using `rho=reshape(full(rho),[size(rho,1)/prod(dims) prod(dims) stack_size])`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(rho,dims)`.
  - Representative operation: `if (~isnumeric(rho))`.
  - Representative operation: `error('rho must be a numeric array.')`.

## Parameters / inputs

- rho -Fokker-Planck state vector
- dims -spatial dimensions of the
- Fokker-Planck problem, a
- vector of positive integers

## Outputs

- rho -Liouville space state vector

## Implementation structure

- Returns the average of the spin state vector across the
- spatial dimensions of the sample. Syntax:
- rho=fpl2rho(rho,dims)
- rho -Fokker-Planck state vector
- dims -spatial dimensions of the
- Fokker-Planck problem, a
- vector of positive integers
- rho -Liouville space state vector
- Check consistency
- Find out the stack size
- Expose the spin dimension (no ND sparse support yet)
- Average over the spatial coordinates

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `squeeze()`, `any()`, `space()`.
