# kernel/overloads/@opium/mtimes.m

- Signature: `c=mtimes(a,b)`

## Purpose

Matrix products involving an OPIUM object. Syntax: c=mtimes(a,b)

## Physical / mathematical content

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
