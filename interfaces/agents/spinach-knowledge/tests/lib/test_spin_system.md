# tests/lib/test_spin_system.m

- Signature: `spin_system=test_spin_system(sys,inter,bas)`

## Purpose

Builds a Spinach spin system and basis with the quiet settings used by regression tests.

## Parameters / inputs

- `sys` - Spinach system specification.
- `inter` - Spinach interaction specification.
- `bas` - Spinach basis specification.

## Outputs

- `spin_system` - Spinach spin-system object with the requested basis.

## Implementation structure

- Sets `sys.output='hush'`, adds `hygiene` to `sys.disable` without duplicates, sets `sys.parallel={'local',1}` and `sys.parprops={}`, then calls `create(sys,inter)` and `basis(spin_system,bas)`.
