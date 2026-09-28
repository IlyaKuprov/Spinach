# kernel/overloads/@polyadic/minus.m

- Signature: `a=minus(a,b)`

## Purpose

Polyadic subtraction operation. Does not perform the subtraction immediately, but instead stores the operands as a sum of unopened Kronecker products. Syntax: c=minus(a,b)

## Physical / mathematical content

Subtraction is represented by negating the second operand and adding it to the first; the polyadic representation avoids immediately expanding the Kronecker products.

## Numerical / algorithmic content

## Parameters / inputs

- `a`, `b`: polyadic objects

## Outputs

- `c`: polyadic object
- Note: use this operation sparingly—the subtractions are buffered, and all subsequent operations will be slower.

## Implementation structure

- Polyadic subtraction operation. Does not perform the actual sub-
- traction, but instead stores the operands as a sum of unopened
- Kronecker products. Syntax:
- c=minus(a,b)
- a,b -polyadic objects
- c -polyadic object
- Note: use this operation sparingly -the subtractions are simply
- buffered, and all subsequent operations will be slower.
- Just call plus
- The 1958 Fourier transform NMR article by Morozov, Melnikov and
- Skripov only came to light during a patent dispute between Bruker
- and Varian. The Nobel Prize winning paper by Ernst and Anderson
