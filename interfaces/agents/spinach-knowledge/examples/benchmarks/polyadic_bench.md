# examples/benchmarks/polyadic_bench.m

- MATLAB implementation: [examples/benchmarks/polyadic_bench.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/polyadic_bench.m)

polyadic_bench() is a computational microbenchmark, not a spin-system simulation. It compares vector multiplication by a three-factor Spinach polyadic with multiplication by an expanded matrix, separately for full and sparse random complex factors. Run with MATLAB and Spinach's polyadic implementation on the path; it has no external data-file input.

The script uses 100 samples (n_stats=100), three matrices per polyadic, and a maximum randomly selected square-factor dimension of 32. In the full case each factor has a dimension drawn independently from 1 through 32 and entries generated from independent real and imaginary Gaussian draws. A random complex vector of matching length is multiplied by P, then by full(P). In the sparse case, each factor is generated with sprandn at requested density 5/dim, multiplied by the phase exp(1i/3); a random complex vector is tested with P and then with inflate(P).

For each of the four cases, elapsed time is converted to milliseconds. The script prints the sample mean and std(samples)/sqrt(100) (labelled “stdev” in its output) as the uncertainty estimate. It does not set a random seed, save the timing arrays, make plots, or compare multiplication outputs for equality.

A source-specific distinction matters when reusing this benchmark: the full representation is formed with full(P), whereas the sparse representation is formed with inflate(P). These are the two different expansion paths whose multiplication costs it measures. Reported times depend on the MATLAB/Spinach build and host; the source supplies no expected benchmark values.
