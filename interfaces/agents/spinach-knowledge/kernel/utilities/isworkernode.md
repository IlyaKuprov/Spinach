# kernel/utilities/isworkernode.m

- Signature: `answer=isworkernode()`

## Purpose

Reports whether the current execution is on a parallel-pool worker. Spinach uses this query to select serial versions of algorithms when already running inside a parallel loop.

## Parameters / inputs

None.

## Outputs

- `answer` - true if running on a parallel worker process.

## Implementation

The function returns the result of MATLAB's undocumented `parallel.internal.pool.isPoolWorker()` function.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=isworkernode.m)
