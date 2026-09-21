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
