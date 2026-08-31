# kernel/overloads/@polyadic/allfinite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/allfinite.m`
- Signature: `answ=allfinite(p)`
- Total lines: 51

## Purpose

Returns true if none of the elements of the polyadic are Inf or NaN. Syntax: answ=allfinite(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check the core array; implemented by `for n=1:numel(p.cores)`.
- Lines 30-31: Check prefix and suffix arrays; implemented by `for n=1:numel(p.prefix)`.
- Lines 42-43: All finite; implemented by `answ=true()`.

### Control flow inferred from the code

- Line 22: `for` loop over `n=1:numel(p.cores)`.
- Line 23: `for` loop over `k=1:numel(p.cores{n})`.
- Line 24: conditional branch on `~allfinite(p.cores{n}{k})`.
- Line 31: `for` loop over `n=1:numel(p.prefix)`.
- Line 32: conditional branch on `~allfinite(p.prefix{n})`.
- Line 36: `for` loop over `n=1:numel(p.suffix)`.
- Line 37: conditional branch on `~allfinite(p.suffix{n})`.

### Key state/data transformations

- Lines 25: computes `answ` using `answ=false; return`.

## Parameters / inputs

- p -a polyadic object

## Outputs

- answ -logical true if all numeric data
- in the polyadic object is finite

## Implementation structure

- Returns true if none of the elements of the polyadic are Inf
- or NaN. Syntax:
- answ=allfinite(p)
- p -a polyadic object
- answ -logical true if all numeric data
- in the polyadic object is finite
- Check the core array
- Check prefix and suffix arrays
- All finite
- Beauty is the first test: there is no permanent
- place in the world for ugly mathematics.
- G.H. Hardy

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `true()`.
