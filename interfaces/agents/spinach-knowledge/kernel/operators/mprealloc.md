# kernel/operators/mprealloc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/mprealloc.m`
- Signature: `A=mprealloc(spin_system,nnzpc)`
- Total lines: 72

## Purpose

Preallocates an operator in the current basis. Syntax: A=mprealloc(spin_system,nnzpc)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(spin_system,nnzpc)`.
- Lines 24-25: Do the math; implemented by `switch spin_system.bas.formalism`.
- Lines 29-30: Create a zero Liouville space matrix operator; implemented by `problem_dim=size(spin_system.bas.basis,1)`.
- Lines 35-36: Create a zero Hilbert space matrix operator; implemented by `problem_dim=prod(spin_system.comp.mults)`.
- Lines 41-42: Create a zero Liouville space matrix operator; implemented by `problem_dim=prod(spin_system.comp.mults.^2)`.
- Lines 47-48: Complain and bomb out; implemented by `error('unknown formalism specification.')`.

### Control flow inferred from the code

- Line 25: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `{'zeeman-wavef','zeeman-hilb'}`, `'zeeman-liouv'`.

### Key state/data transformations

- Lines 30: computes `problem_dim` using `problem_dim=size(spin_system.bas.basis,1)`.
- Lines 31: computes `A` using `A=spalloc(problem_dim,problem_dim,nnzpc*problem_dim)`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(spin_system,nnzpc)`. In 2006, Oxford's Magdalen College, where Erwin Schrodinger was a Fellow between 1933 and 1936, received a sum of money from a benefactor towards
  - Representative operation: `if (~isnumeric(nnzpc))||(~isreal(nnzpc))||(~isscalar(nnzpc))||(mod(nnzpc,1)~=0)`.
  - Representative operation: `error('nnzpc parameter must be a positive real integer.')`.

## Parameters / inputs

- nnzpc -expected number of non-zeros per column

## Outputs

- A -all-zero sparse matrix of the appropriate
- dimension with room for the specified num-
- ber of non-zeroes

## Implementation structure

- Preallocates an operator in the current basis. Syntax:
- A=mprealloc(spin_system,nnzpc)
- nnzpc - expected number of non-zeros per column
- A - all-zero sparse matrix of the appropriate
- dimension with room for the specified num-
- ber of non-zeroes
- Check consistency
- Do the math
- Create a zero Liouville space matrix operator
- Create a zero Hilbert space matrix operator
- Complain and bomb out
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `isscalar()`, `isfield()`.
