# kernel/utilities/poolsize.m

## Purpose

Returns the current parallel pool size in Spinach ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/poolsize.m)).

## Behaviour

- Syntax: `n=poolsize()`.
- The function obtains the current parallel pool handle with `gcp('nocreate')`.
- If no pool exists, `n` is set to `0`.
- Otherwise, `n` is set to `p.NumWorkers`, the number of workers in the current parallel pool.
- When invoked from inside `parfor`, `spmd`, or an asynchronous parallel job, the function returns zero.

## Inputs and outputs

**Inputs**

- None.

**Outputs**

- `n` — number of workers in the current parallel pool.

## References

- [Spinach Wiki: poolsize.m](https://spindynamics.org/wiki/index.php?title=poolsize.m)
- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/poolsize.m)
