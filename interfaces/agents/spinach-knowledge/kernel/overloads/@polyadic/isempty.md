# kernel/overloads/@polyadic/isempty.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/isempty.m`
- Signature: `answer=isempty(p)`
- Total lines: 33

## Purpose

Returns true for polyadics that represent a matrix with a zero dimension. Syntax: answer=isempty(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -a logical value

## Implementation structure

- Returns true for polyadics that represent a matrix with a
- zero dimension. Syntax:
- answer=isempty(p)
- p -a polyadic object
- answer -a logical value
- Check the size
- Q: Why do sumo wrestlers shave their legs and armpits?
- A: To make sure people can tell them apart from feminists.
- A "festive season" cracker

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `any()`, `true()`, `false()`.
