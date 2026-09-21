# kernel/overloads/@opium/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@opium/mtimes.m`
- Signature: `c=mtimes(a,b)`
- Total lines: 85

## Purpose

Matrix products involving an OPIUM object. Syntax: c=mtimes(a,b)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Parameters / inputs

- a,b -opia or numerical arrays

## Outputs

- c -multiplication result

## Implementation structure

- Matrix products involving an OPIUM object. Syntax:
- c=mtimes(a,b)
- a,b -opia or numerical arrays
- c -multiplication result
- When A is a scalar
- Return opium multiplied by A
- When A is not a scalar
- Check dimension
- Return A multiplied by opium
- When B is a scalar
- When B is not a scalar
- When both are opia

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`.
