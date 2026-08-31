# kernel/operators/sin_tran.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/sin_tran.m`
- Signature: `A=sin_tran(dim)`
- Total lines: 67

## Purpose

Single transition operators, spanning the space of matri- ces of the specified dimension. The set is returned as a cell array of sparse matrices using serpentine indexing where the position in the cell array maps in the follow- ing way to the location of the single non-zero: (1) (3) (6) (10) (2) (5) (9) (13) (4) (8) (12) (15) (7) (11) (14) (16) and likewise for larger matrices. Syntax: A=sin_tran(dim)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(dim)`.
- Lines 36-37: Empty array; implemented by `A=cell(dim^2,1)`.
- Lines 39-40: Fill the array; implemented by `parfor n=1:dim^2`.
- Lines 42-43: Unit element indices; implemented by `[k,q]=lin2kq(dim,n,1)`.
- Lines 45-46: Matrix construction; implemented by `A{n}=sparse(k,q,1,dim,dim)`.
- Lines 48-49: Complex type; implemented by `A{n}=complex(A{n})`.

### Control flow inferred from the code

- Line 40: `parfor` loop over `n=1:dim^2`.

### Key state/data transformations

- Lines 37: computes `A` using `A=cell(dim^2,1)`.
- Lines 43: computes `[k,q]` using `[k,q]=lin2kq(dim,n,1)`.
- Lines 46: computes `A{n}` using `A{n}=sparse(k,q,1,dim,dim)`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(dim)`. Have nothing in your house that you do not know to be useful or believe to be beautiful.
  - Representative operation: `if (~isnumeric(dim))||(~isscalar(dim))|| (~isreal(dim))||(dim<1)||(mod(dim,1)~=0)`.
  - Representative operation: `(~isreal(dim))||(dim<1)||(mod(dim,1)~=0)`.

## Parameters / inputs

- dim -dimension of the matrices

## Outputs

- A -a cell array of matrices, structured
- as described above; matrices are re-
- turned as complex to avoid expensive
- reallocations later

## Implementation structure

- Single transition operators, spanning the space of matri-
- ces of the specified dimension. The set is returned as a
- cell array of sparse matrices using serpentine indexing
- where the position in the cell array maps in the follow-
- ing way to the location of the single non-zero:
- (1) (3) (6) (10)
- (2) (5) (9) (13)
- (4) (8) (12) (15)
- (7) (11) (14) (16)
- and likewise for larger matrices. Syntax:
- A=sin_tran(dim)
- dim -dimension of the matrices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lin2kq()`, `complex()`, `isscalar()`.
