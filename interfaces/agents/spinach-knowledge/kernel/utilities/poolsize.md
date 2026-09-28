# kernel/utilities/poolsize.m

- Signature: `n=poolsize()`

## Purpose

Returns the number of workers in the current parallel pool, or `0` when no pool exists.

## Parameters / inputs

None.

## Output

- `n` — number of workers in the current parallel pool.

## Implementation structure

The function calls `gcp('nocreate')` to query an existing pool without creating one. If the result is empty it returns `0`; otherwise it returns the pool's `NumWorkers`. The source documents that calls from `parfor`, `spmd`, or an asynchronous parallel job return `0`.

## Reference

- <https://spindynamics.org/wiki/index.php?title=poolsize.m>
