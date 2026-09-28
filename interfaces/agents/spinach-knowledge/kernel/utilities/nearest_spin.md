# kernel/utilities/nearest_spin.m

- Signature: `[k,d]=nearest_spin(spin_system,n)`

## Purpose

Returns the index of the nearest spin to the one specified. Only spins for which Cartesian coordinates are available are considered. Syntax: k=nearest_spin(spin_system,n)

## Physical / mathematical content

- Finds the closest other spin to spin `n` using the spins’ Cartesian coordinates.

## Numerical / algorithmic content

- Scans other spins with available coordinates and compares Euclidean 2-norm distances; it excludes `n` and errors if no other spin has coordinates.

## Parameters / inputs

- n -index of the spin in question

## Outputs

- k -index of the nearest spin
- d -distance to the nearest spin, Angstrom

## Implementation structure

- Validates that `n` is an existing positive-integer spin index with coordinates, then retains the smallest distance among eligible spins and returns its index and distance.
