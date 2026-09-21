# kernel/overloads/@polyadic/allfinite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/allfinite.m`
- Signature: `answ=allfinite(p)`
- Total lines: 51

## Purpose

Returns true if none of the elements of the polyadic are Inf or NaN. Syntax: answ=allfinite(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

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
