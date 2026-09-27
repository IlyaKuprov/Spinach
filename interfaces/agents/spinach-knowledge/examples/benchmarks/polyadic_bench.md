# examples/benchmarks/polyadic_bench.m

- Signature: `polyadic_bench()`

## Purpose

A benchmark for the polyadic object.

## Numerical / algorithmic content

- Compares multiplication by three-factor polyadics with multiplication by their full or inflated matrix forms, for both full and sparse random complex factors. The benchmark reports mean runtimes and standard errors over 100 samples.

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
