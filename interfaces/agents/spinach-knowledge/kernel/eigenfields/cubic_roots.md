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
