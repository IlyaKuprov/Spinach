# kernel/conventions/transforms/hartree2joule.m

- Signature: `energy=hartree2joule(energy)`

## Purpose

Converts Hartree energy units into J/mol. A Hartree is twice the ground state ionisation energy of the hydrogen atom. Syntax: energy=hartree2joule(energy)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- energy -a numerical array of energies in
- Hartree units

## Outputs

- energy -a numerical array of energies in
- Joules

## Implementation structure

- Converts Hartree energy units into J/mol. A Hartree is twice the
- ground state ionisation energy of the hydrogen atom. Syntax:
- energy=hartree2joule(energy)
- energy -a numerical array of energies in
- Hartree units
- Joules
- Check consistency
- Perform the conversion
- Consistency enforcement
- "We only need to be lucky once. You need to be
- lucky every time."
- The IRA to Margaret Thatcher, after
