# kernel/operators/twospinist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/twospinist.m`
- Signature: `T=twospinist(spin_system,spin_a,spin_b,indices,type)`
- Total lines: 148

## Purpose

Two-spin irreducible spherical tensor operators. Syntax: T=twospinist(spin_system,spin_a,spin_b,indices,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,spin_a,spin_b,indices,type)`.
- Lines 41-42: Indexing switches; implemented by `switch indices(1)`.
- Lines 46-47: Three first rank operators; implemented by `switch indices(2)`.
- Lines 66-67: Complain and bomb out; implemented by `error('incorrect spherical tensor indices.')`.
- Lines 73-74: Five second rank operators; implemented by `switch indices(2)`.

### Control flow inferred from the code

- Line 42: dispatches on `indices(1)`; cases `1`, `+1`, `0`, `-1`.
- Line 47: dispatches on `indices(2)`; cases `+1`, `0`, `-1`.
- Line 74: dispatches on `indices(2)`; cases `+2`, `+1`, `0`, `-1`, `-2`.

### Key state/data transformations

- Lines 51-52: computes `T` using `T=-(1/2)*(operator(spin_system,{'L+','Lz'},{spin_a spin_b},type,'csc')- operator(spin_system,{'Lz','L+'},{spin_a spin_b},type,'csc'))`.

### Local helper functions

- Line 117: `grumble()` — `function grumble(spin_system,spin_a,spin_b,indices,type)`.
  - Representative operation: `if (~isnumeric(spin_a))||(~isscalar(spin_a))||(~isreal(spin_a))|| (spin_a<1)||(spin_a>spin_system.comp.nspins)||(mod(spin_a,1)~=0)`.
  - Representative operation: `(spin_a<1)||(spin_a>spin_system.comp.nspins)||(mod(spin_a,1)~=0)`.

## Parameters / inputs

- spin_a -number of the first spin
- spin_b -number of the second spin
- indices -two-element vector [L,M] containing
- rank and projection index; L=1,2 are
- available
- In Liouville space, type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator
- In Hilbert space calculations, the type parameter is ig-
- nored, and the operator itself is always returned.

## Outputs

- T -irreducible spherical tensor operator

## Implementation structure

- Two-spin irreducible spherical tensor operators. Syntax:
- T=twospinist(spin_system,spin_a,spin_b,indices,type)
- spin_a -number of the first spin
- spin_b -number of the second spin
- indices -two-element vector [L,M] containing
- rank and projection index; L=1,2 are
- available
- In Liouville space, type can be set to:
- 'left' -produces left side product superoperator
- 'right' -produces right side product superoperator
- 'comm' -produces commutation superoperator (default)
- 'acomm' -produces anticommutation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `indices()`, `operator()`, `isscalar()`, `ischar()`, `ismember()`.
