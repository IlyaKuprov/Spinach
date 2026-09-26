# etc/diamond_defects/diamond_ov0.m

- Signature: `[sys,inter]=diamond_ov0(parameters)`

## Purpose

Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_ov0(parameters) Magnetic parameters from Tables 9-2 and 9-3 of: B.L. Cann, Magnetic Resonance Studies of Point Defects in Diamond, PhD thesis, University of Warwick (2009), The zero-field splitting is confirmed at 4 K, and the defect is assigned to OV0, in: S. Mukherjee et al., Phys. Rev. B 114, 074105 (2026), 

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- a structure (parameters.*) with the following field:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system
- for diamond. Syntax:
- [sys,inter]=diamond_ov0(parameters)
- Magnetic parameters from Tables 9-2 and 9-3 of:
- B.L. Cann, Magnetic Resonance Studies of Point Defects in
- Diamond, PhD thesis, University of Warwick (2009),
- The zero-field splitting is confirmed at 4 K, and the defect is
- assigned to OV0, in:
- S. Mukherjee et al., Phys. Rev. B 114, 074105 (2026),
- The centre has S=1 and C3v symmetry; the electron Zeeman and zero-
- field splitting tensors are axial about the trigonal axis, and the
- rhombicity is zero within the experimental error. Oxygen is left out
