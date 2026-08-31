# kernel/overloads/@polyadic/minus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/minus.m`
- Signature: `a=minus(a,b)`
- Total lines: 38

## Purpose

Polyadic subtraction operation. Does not perform the actual sub- traction, but instead stores the operands as a sum of unopened Kronecker products. Syntax: c=minus(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Just call plus; implemented by `a=plus(a,(-1)*b)`.

### Key state/data transformations

- Lines 25: computes `a` using `a=plus(a,(-1)*b)`.

## Parameters / inputs

- a,b -polyadic objects

## Outputs

- c -polyadic object
- Note: use this operation sparingly -the subtractions are simply
- buffered, and all subsequent operations will be slower.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `plus()`.
