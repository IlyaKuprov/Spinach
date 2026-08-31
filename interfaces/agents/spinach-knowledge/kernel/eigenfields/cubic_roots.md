# kernel/eigenfields/cubic_roots.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/cubic_roots.m`
- Signature: `root_list=cubic_roots(poly_coeffs,root_tol)`
- Total lines: 86

## Purpose

Real roots of a cubic polynomial in the unit interval. Syntax: root_list=cubic_roots(poly_coeffs,root_tol)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(poly_coeffs,root_tol)`.
- Lines 25-26: Normalise polynomial coefficients; implemented by `poly_coeffs=poly_coeffs(:).'`.
- Lines 34-35: Drop leading numerical zeros; implemented by `lead_idx=find(abs(poly_coeffs)>root_tol,1,'first')`.
- Lines 42-43: Find real roots inside the unit interval; implemented by `root_list=roots(poly_coeffs)`.
- Lines 51-52: Merge numerically coincident roots; implemented by `if ~isempty(root_list)`.

### Control flow inferred from the code

- Line 28: conditional branch on `poly_scale==0`.
- Line 36: conditional branch on `isempty(lead_idx)`.
- Line 52: conditional branch on `~isempty(root_list)`.

### Key state/data transformations

- Lines 26: computes `poly_coeffs` using `poly_coeffs=poly_coeffs(:).'`.
- Lines 27: computes `poly_scale` using `poly_scale=max(abs(poly_coeffs))`.
- Lines 29: computes `root_list` using `root_list=[]`.
- Lines 35: computes `lead_idx` using `lead_idx=find(abs(poly_coeffs)>root_tol,1,'first')`.

### Local helper functions

- Line 59: `grumble()` — `function grumble(poly_coeffs,root_tol)`. You can always tell when a public figure has said something with the ring
  - Representative operation: `if (~isnumeric(poly_coeffs))||(~isreal(poly_coeffs))|| (numel(poly_coeffs)~=4)||any(~isfinite(poly_coeffs(:)))`.
  - Representative operation: `(numel(poly_coeffs)~=4)||any(~isfinite(poly_coeffs(:)))`.

## Parameters / inputs

- poly_coeffs -four real coefficients [a b c d] of
- a*x^3+b*x^2+c*x+d
- root_tol -positive real root filtering tolerance

## Outputs

- root_list -sorted row vector of real roots in [0,1]

## Implementation structure

- Real roots of a cubic polynomial in the unit interval. Syntax:
- root_list=cubic_roots(poly_coeffs,root_tol)
- poly_coeffs -four real coefficients [a b c d] of
- a*x^3+b*x^2+c*x+d
- root_tol -positive real root filtering tolerance
- root_list -sorted row vector of real roots in [0,1]
- Check consistency
- Normalise polynomial coefficients
- Drop leading numerical zeros
- Find real roots inside the unit interval
- Merge numerically coincident roots
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `poly_coeffs()`, `roots()`, `root_list()`, `diff()`, `any()`, `isscalar()`.
