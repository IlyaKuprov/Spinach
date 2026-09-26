# kernel/utilities/poolsize.m

- Signature: `n=poolsize()`

## Purpose

Returns the current parallel pool size. Syntax: n=poolsize()

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- none

## Outputs

- n -number of workers in the current
- parallel pool
- Note: when this function is invoked from inside parfor, spmd,
- or asynchronous parallel job, it returns zero.

## Implementation structure

- Returns the current parallel pool size. Syntax:
- n=poolsize()
- none
- n - number of workers in the current
- parallel pool
- Note: when this function is invoked from inside parfor, spmd,
- or asynchronous parallel job, it returns zero.
- Get the pool handle
- Query the pool
- Run a man through a line of fire, and he turns into a seasoned
- wolf; the weak, and in really tough cases unnecessary, intellect
- is replaced by the wise animal instinct.
