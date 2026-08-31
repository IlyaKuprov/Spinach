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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 70-71: Check consistency; implemented by `grumble(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`.
- Lines 73-76: Get traingle subdivision midpoints; implemented by `[r12,r23,r31]=sphtrsubd(triangle(1).xyz, triangle(2).xyz, triangle(3).xyz)`.
- Lines 78-79: First subdivision midpoint eigenset; implemented by `[phi,elev]=cart2sph(r12(1),r12(2),r12(3))`.
- Lines 87-88: Second subdivision midpoint eigenset; implemented by `[phi,elev]=cart2sph(r23(1),r23(2),r23(3))`.
- Lines 96-97: Third subdivision midpoint eigenset; implemented by `[phi,elev]=cart2sph(r31(1),r31(2),r31(3))`.
- Lines 105-106: Eigensets for the vertices (three each) of the four subdivision traingles; implemented by `triangle_a(1)=triangle(1); triangle_a(2)=eigenset12; triangle_a(3)=eigenset31`.
- Lines 111-115: Compute the subdivided integral; implemented by `spec_sub=trint(parameters,triangle_a)+ trint(parameters,triangle_b)+ trint(parameters,triangle_c)+ trint(parameters,triangle_d)`.
- Lines 117-118: Compute the direct integral; implemented by `spec_dir=trint(parameters,triangle)`.
- Lines 120-121: Simpson-Richardson correction; implemented by `spec_sim=(4*spec_sub-spec_dir)/3`.
- Lines 123-124: If the accuracy is insufficient, recurse; implemented by `if norm(spec_dir-spec_sub,2)>parameters.int_tol`.
- Lines 126-127: Compute the four triangles of the subdivision via asynchronous recursion; implemented by `spec_a=parfeval(@voitlander,1,spin_system,parameters,triangle_a,Ic,Iz,Qc,Qz,Hmw)`.
- Lines 132-134: Retrieve the work; implemented by `spec=fetchOutputs(spec_a)+fetchOutputs(spec_b)+ fetchOutputs(spec_c)+fetchOutputs(spec_d)`.
- Lines 138-139: Good enough; implemented by `spec=spec_sim`.

### Control flow inferred from the code

- Line 124: conditional branch on `norm(spec_dir-spec_sub,2)>parameters.int_tol`.

### Key state/data transformations

- Lines 74-76: computes `[r12,r23,r31]` using `[r12,r23,r31]=sphtrsubd(triangle(1).xyz, triangle(2).xyz, triangle(3).xyz)`.
- Lines 79: computes `[phi,elev]` using `[phi,elev]=cart2sph(r12(1),r12(2),r12(3))`.
- Lines 80: computes `parameters12` using `parameters12=parameters`.
- Lines 81: computes `parameters12.orientation` using `parameters12.orientation=[0, pi/2-elev, phi]`.
- Lines 82: computes `Hz` using `Hz=Iz+orientation(Qz,parameters12.orientation)`.
- Lines 83: computes `Hc` using `Hc=Ic+orientation(Qc,parameters12.orientation)`.
- Lines 84: computes `eigenset12` using `eigenset12=eigenfields(spin_system,parameters12,Hz,Hc,Hmw)`.
- Lines 85: computes `eigenset12.xyz` using `eigenset12.xyz=[r12(1); r12(2); r12(3)]`.
- Lines 89: computes `parameters23` using `parameters23=parameters`.
- Lines 90: computes `parameters23.orientation` using `parameters23.orientation=[0, pi/2-elev, phi]`.
- Lines 93: computes `eigenset23` using `eigenset23=eigenfields(spin_system,parameters23,Hz,Hc,Hmw)`.
- Lines 94: computes `eigenset23.xyz` using `eigenset23.xyz=[r23(1); r23(2); r23(3)]`.
- Lines 98: computes `parameters31` using `parameters31=parameters`.
- Lines 99: computes `parameters31.orientation` using `parameters31.orientation=[0, pi/2-elev, phi]`.
- Lines 102: computes `eigenset31` using `eigenset31=eigenfields(spin_system,parameters31,Hz,Hc,Hmw)`.
- Lines 103: computes `eigenset31.xyz` using `eigenset31.xyz=[r31(1); r31(2); r31(3)]`.
- Lines 106: computes `triangle_a(1)` using `triangle_a(1)=triangle(1); triangle_a(2)=eigenset12; triangle_a(3)=eigenset31`.
- Lines 107: computes `triangle_b(1)` using `triangle_b(1)=eigenset12; triangle_b(2)=triangle(2); triangle_b(3)=eigenset23`.

### Local helper functions

- Line 146: `trint()` — `function spec=trint(parameters,triangle)`. Preallocate the spectrum
  - Representative operation: `spec=zeros(size(parameters.b_axis),'like',1i)`.
  - Representative operation: `if isempty(triangle(1).ti)|| isempty(triangle(2).ti)|| isempty(triangle(3).ti)`.
- Line 282: `grumble()` — `function grumble(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)`.
  - Representative operation: `if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))|| (~isfield(spin_system.bas,'formalism'))`.
  - Representative operation: `(~isfield(spin_system.bas,'formalism'))`.

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
