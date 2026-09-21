# kernel/conventions/transforms/ham2nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ham2nqi.m`
- Signature: `[omega,Q]=ham2nqi(H)`
- Total lines: 94

## Purpose

Converts a single-spin Hamiltonian back into the Zeeman and quadrupolar interaction parameters that had been used to generate it. Syntax: [omega,Q]=ham2nqi(H)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H -single-spin Hamiltonian written in
- the Zeeman basis for a spin of any
- multiplicity

## Outputs

- omega -Larmor frequencies, rad/s
- Q -symmetric traceless quadrupolar
- coupling tensor, rad/s
- The outputs are returned such that:
- H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz +
- + [Sx Sy Sz]*Q*[Sx Sy Sz].';
- An error is produced if the Hamilonian contains
- any terms (for example, cubic) beyond those, or
- if it is not Hermitian and traceless.

## Implementation structure

- Converts a single-spin Hamiltonian back into the
- Zeeman and quadrupolar interaction parameters that
- had been used to generate it. Syntax:
- [omega,Q]=ham2nqi(H)
- H -single-spin Hamiltonian written in
- the Zeeman basis for a spin of any
- multiplicity
- omega -Larmor frequencies, rad/s
- Q -symmetric traceless quadrupolar
- coupling tensor, rad/s
- The outputs are returned such that:
- H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz +

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `irr_sph_ten()`, `pauli()`, `omega()`, `rank2()`, `sphten2mat()`, `ishermitian()`, `eps()`.
