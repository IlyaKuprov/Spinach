# kernel/overloads/@polyadic/polyadic.m

- Signature: `p=polyadic(cores)`

## Purpose

Creates an object of a polyadic class. Syntax: p=polyadic(cores) A polyadic is a matrix formed by a Kronecker product, with those products stored unopened. For example, cores={{A,B,C},{D,E}} corresponds to A(x)B(x)C + D(x)E matrix. Any multiplicative acti- on by this matrix may be computed without opening the Kronecker products. This can save orders of magnitude in CPU time.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- cores -a cell array of cell arrays of matrices
- whose Konecker products make up the mat-
- rix of interest.

## Outputs

- p -a polyadic representation of a matrix
- that behaves in many respects like the
- matrix it represents.
- Note: nested polyadics are permitted -the input matrices may be
- polyadics themselves.

## Implementation structure

- Creates an object of a polyadic class. Syntax:
- p=polyadic(cores)
- A polyadic is a matrix formed by a Kronecker product, with those
- products stored unopened. For example,
- cores={{A,B,C},{D,E}}
- corresponds to A(x)B(x)C + D(x)E matrix. Any multiplicative acti-
- on by this matrix may be computed without opening the Kronecker
- products. This can save orders of magnitude in CPU time.
- cores - a cell array of cell arrays of matrices
- whose Konecker products make up the mat-
- rix of interest.
- p - a polyadic representation of a matrix
