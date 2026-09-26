# kernel/includes/parallel_profiler_start.m

- Signature: `(script file)`

## Purpose

An include that starts profiling infrastructure around parallel stages. Should be invoked just before a parfor or an spmd.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- An include that starts profiling infrastructure around parallel
- stages. Should be invoked just before a parfor or an spmd.
- Brief parallel profiler start
- Detailed parallel profiler start
- In late 1700s, a teacher in a German school asked a kid to
- sum up the numbers from 1 to 100 as a punishment for misbe-
- having. The teacher was astonished when the kid solved the
- problem in seconds:
- S = 1 + 2 + ... + 100
- S = 100 + 99 + ... + 1
- --------------------------
- 2S = 101 + 101 + ... + 101 => S = 101*100/2 = 5050
