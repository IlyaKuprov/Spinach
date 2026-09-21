# kernel/eigenfields/rootmatch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/rootmatch.m`
- Signature: `[idx1,idx2,idx3]=rootmatch(field1,field2,field3,...`
- Total lines: 214

## Purpose

Global order-preserving root matching between three magnetic field root lists. The routine returns indices into the input lists that identify matching roots with minimum field-continuation cost. Syntax: [idx1,idx2,idx3]=rootmatch(field1,field2,field3,... edge12,edge23,edge31)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- field1 -real vector of roots at the first triangle vertex
- field2 -real vector of roots at the second triangle vertex
- field3 -real vector of roots at the third triangle vertex
- edge12 -positive distance between the first and the second
- triangle vertices
- edge23 -positive distance between the second and the third
- triangle vertices
- edge31 -positive distance between the third and the first
- triangle vertices

## Outputs

- idx1 -indices of matched roots in field1
- idx2 -indices of matched roots in field2
- idx3 -indices of matched roots in field3

## Implementation structure

- Global order-preserving root matching between three magnetic field root
- lists. The routine returns indices into the input lists that identify
- matching roots with minimum field-continuation cost. Syntax:
- [idx1,idx2,idx3]=rootmatch(field1,field2,field3,...
- edge12,edge23,edge31)
- field1 -real vector of roots at the first triangle vertex
- field2 -real vector of roots at the second triangle vertex
- field3 -real vector of roots at the third triangle vertex
- edge12 -positive distance between the first and the second
- triangle vertices
- edge23 -positive distance between the second and the third
- edge31 -positive distance between the third and the first

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `field1()`, `field2()`, `field3()`, `inf()`, `match_cost()`, `uint8()`, `match_count()`, `prev_move()`, `idx1()`, `ord1()`, `idx2()`, `ord2()`, `idx3()`, `ord3()`, `fliplr()`.
