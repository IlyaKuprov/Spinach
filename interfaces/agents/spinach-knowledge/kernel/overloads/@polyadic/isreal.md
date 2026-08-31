# kernel/overloads/@polyadic/isreal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/isreal.m`
- Signature: `answer=isreal(p)`
- Total lines: 48

## Purpose

Returns true if the polyadic representation is real. Syntax: answer=isreal(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Check the core array; implemented by `for n=1:numel(p.cores)`.
- Lines 28-29: Check prefix and suffix arrays; implemented by `for n=1:numel(p.prefix)`.
- Lines 40-41: All data is real; implemented by `answer=true()`.

### Control flow inferred from the code

- Line 20: `for` loop over `n=1:numel(p.cores)`.
- Line 21: `for` loop over `k=1:numel(p.cores{n})`.
- Line 22: conditional branch on `~isreal(p.cores{n}{k})`.
- Line 29: `for` loop over `n=1:numel(p.prefix)`.
- Line 30: conditional branch on `~isreal(p.prefix{n})`.
- Line 34: `for` loop over `n=1:numel(p.suffix)`.
- Line 35: conditional branch on `~isreal(p.suffix{n})`.

### Key state/data transformations

- Lines 23: computes `answer` using `answer=false; return`.

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -true if all numeric data in the object is real

## Implementation structure

- Returns true if the polyadic representation is real. Syntax:
- answer=isreal(p)
- p -a polyadic object
- answer -true if all numeric data in the object is real
- Check the core array
- Check prefix and suffix arrays
- All data is real
- A little inaccuracy sometimes saves a ton of explanation.
- H.H. Munro

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `true()`.
