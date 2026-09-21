# examples/benchmarks/polyadic_bench.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/polyadic_bench.m`
- Signature: `polyadic_bench()`
- Total lines: 100

## Purpose

A benchmark for the polyadic object.

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

## Implementation structure

- A benchmark for the polyadic object.
- Statistics parameters
- % Full matrix benchmark
- Result array
- Full matrix statistics loop
- Update the user
- Get random full complex matrices
- Form a polyadic
- Get a random full complex vector
- Time polyadic multiplication
- Inflate the polyadic
- Time flat multiplication

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2str()`, `randi()`, `polyadic()`, `runtimes_full_poly()`, `runtimes_full_flat()`, `sprandn()`, `runtimes_sparse_poly()`, `inflate()`, `runtimes_sparse_flat()`, `std()`.
