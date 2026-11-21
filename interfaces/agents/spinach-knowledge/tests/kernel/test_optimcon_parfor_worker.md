# tests/kernel/test_optimcon_parfor_worker.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_parfor_worker.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_parfor_worker.m)

## Purpose

Regression test for optimal control launched from inside an external `parfor` loop. When `optimcon()` and `fmaxnewton()` run in the body of the loop, a process-pool worker has no visible pool. `optimcon()` then prepares the problem without per-worker drift slices, and `ensemble()` evaluates every ensemble case serially on that worker. The test requires the result to be the same optimisation the client would perform. It also confirms that `ensemble()` still rejects a problem that `optimcon()` prepared for the client pool when that problem reaches a worker.

## Behaviour

- Opens a two-worker process pool if none exists.
- Defines four distinct one-spin Hilbert-space state-transfer problems, Z to X magnetisation under different frequency offsets. Two of them use power-level ensembles with two and three members.
- Optimises each problem on the client with ten L-BFGS iterations of `fmaxnewton(...,@grape_xy,...)` from a deterministic guess.
- Runs the same four optimisations, `optimcon()` included, inside `parfor`. Each must complete and reproduce its client waveform to a relative tolerance of 1e-10.
- Calls `optimcon()` on the client, passes the prepared problem into `parfor`, and requires `grape_xy()` on every worker to raise the `ensemble()` pool-identity error.

## Inputs and outputs

The test takes no input arguments and returns the standard regression-test `result` structure.

## References

- [Spinach GitHub repository — test_optimcon_parfor_worker.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_parfor_worker.m)
- [Spinach GitHub repository (main page)](https://github.com/IlyaKuprov/Spinach)
