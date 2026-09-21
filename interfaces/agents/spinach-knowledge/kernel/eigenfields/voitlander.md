# kernel/eigenfields/voitlander.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/voitlander.m`
- Signature: `spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`
- Total lines: 421

## Purpose

Adaptively recursed Voitlander integrator. Computes an approximation of an integral of field-swept EPR transition over a spherical triang- le. Syntax: spec=voitlander(spin_system,parameters,... triangle,Ic,Iz,Qc,Qz,Hmw)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `trint()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- triangle(1:3).xyz -Cartesian coordinates of the corners of the
- spherical triangle, unit column vectors
- triangle(1:3).tf -transition fields at the corners of the
- spherical triangle, real column vectors, one
- element per transition
- triangle(1:3).tm -transition moments at the corners of the sphe-
- rical triangle, positive column vectors, one
- element per transition
- triangle(1:3).tw -transition widths at the corners of the sphe-
- rical triangle, positive column vectors, one
- element per transition
- triangle(1:3).pd -energy level population differences at the
- triangle corners, real column vectors, one
- element per transition
- triangle(1:3).ti -transition identity arrays at the triangle
- corners, one row per transition
- triangle(1:3).tj -scaled field-sweep Jacobians at the triangle
- corners, real column vectors, one element per
- transition
- Ic -isotropic part of the coupling Hamiltonian,
- a Hermitian matrix (set retention to 'couplings'
- in assume.m and then call hamiltonian.m)
- Qc -irreducible components of the anisotropic part
- of the coupling Hamiltonian a cell array re-
- turned by hamiltonian.m
- Iz -isotropic part of the Zeeman Hamiltonian, a
- Hermitian matrix (set retention to 'zeeman'
- in assume.m and then call hamiltonian.m) nor-
- malised to 1 Tesla
- Qz -irreducible components of the anisotropic part
- of the Zeeman Hamiltonian a cell array retur-
- ned by hamiltonian.m, normalised to 1 Tesla
- Hmw -perturbation operator, a Hermitian matrix
- parameters.b_axis -a vector of magnetic field values, Tesla
- parameters.int_tol -integration accuracy tolerance

## Outputs

- spec -ESR spectrum integral over the triangle, array
- of the same dimension as parameters.b_axis

## Implementation structure

- Adaptively recursed Voitlander integrator. Computes an approximation
- of an integral of field-swept EPR transition over a spherical triang-
- le. Syntax:
- spec=voitlander(spin_system,parameters,...
- triangle,Ic,Iz,Qc,Qz,Hmw)
- triangle(1:3).xyz -Cartesian coordinates of the corners of the
- spherical triangle, unit column vectors
- triangle(1:3).tf -transition fields at the corners of the
- spherical triangle, real column vectors, one
- element per transition
- triangle(1:3).tm -transition moments at the corners of the sphe-
- rical triangle, positive column vectors, one

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sphtrsubd()`, `triangle()`, `cart2sph()`, `r12()`, `orientation()`, `eigenfields()`, `r23()`, `r31()`, `triangle_a()`, `triangle_b()`, `triangle_c()`, `triangle_d()`, `trint()`, `parfeval()`, `fetchOutputs()`.
