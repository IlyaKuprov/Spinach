# kernel/utilities/isworkernode.m

## Purpose

Returns `true` if the code is currently executing inside a `parfor` or `spmd` block, i.e. on a MATLAB parallel worker process. The function is used in internal Spinach kernel decision making: certain algorithms are switched to their serial versions when the calculation is already running inside a parallel loop.

## Behaviour

- The function takes no arguments and returns a single logical value.
- It calls the undocumented MATLAB internal function `parallel.internal.pool.isPoolWorker()`, which reports whether the current process is a parallel pool worker.

## Inputs and outputs

**Inputs**

- None.

**Outputs**

- `answer` — `true` if running on a parallel worker process, `false` otherwise.

## References

- Source: [kernel/utilities/isworkernode.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isworkernode.m)
- Spinach Wiki: [isworkernode.m](https://spindynamics.org/wiki/index.php?title=isworkernode.m)
