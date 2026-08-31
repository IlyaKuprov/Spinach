# examples/benchmarks/polyadic_bench.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/polyadic_bench.m`
- Signature: `polyadic_bench()`
- Total lines: 100

## Purpose

A benchmark for the polyadic object.

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Statistics parameters; implemented by `n_mats=3`.
- Lines 13-16: % Full matrix benchmark; implemented by `runtimes_full_poly=zeros(n_stats,1)`.
- Lines 15-16: Result array; implemented by `runtimes_full_poly=zeros(n_stats,1)`.
- Lines 19-20: Full matrix statistics loop; implemented by `for n=1:n_stats`.
- Lines 22-23: Update the user; implemented by `disp(['stat loop instance ' num2str(n) '/' num2str(n_stats) ' '])`.
- Lines 25-26: Get random full complex matrices; implemented by `terms=cell(1,n_mats)`.
- Lines 33-34: Form a polyadic; implemented by `P=polyadic({terms})`.
- Lines 36-38: Get a random full complex vector; implemented by `v=1i*randn(size(P,2),1)+ randn(size(P,2),1)`.
- Lines 40-41: Time polyadic multiplication; implemented by `tic; P*v; runtimes_full_poly(n)=1000*toc`.
- Lines 43-44: Inflate the polyadic; implemented by `P=full(P)`.
- Lines 46-47: Time flat multiplication; implemented by `tic; P*v; runtimes_full_flat(n)=1000*toc`.
- Lines 51-54: % Sparse matrix benchmark; implemented by `runtimes_sparse_poly=zeros(n_stats,1)`.
- Lines 53-54: Result array; implemented by `runtimes_sparse_poly=zeros(n_stats,1)`.
- Lines 57-58: Sparse matrix statistics loop; implemented by `for n=1:n_stats`.
- Lines 63-64: Get random sparse complex matrices; implemented by `terms=cell(1,n_mats)`.
- Lines 77-78: Time polyadic multiplication; implemented by `tic; P*v; runtimes_sparse_poly(n)=1000*toc`.
- Lines 80-81: Inflate the polyadic; implemented by `P=inflate(P)`.
- Lines 83-84: Time flat multiplication; implemented by `tic; P*v; runtimes_sparse_flat(n)=1000*toc`.

### Control flow inferred from the code

- Line 20: `for` loop over `n=1:n_stats`.
- Line 27: `for` loop over `k=1:n_mats`.
- Line 58: `for` loop over `n=1:n_stats`.
- Line 65: `for` loop over `k=1:n_mats`.

### Key state/data transformations

- Lines 8: computes `n_mats` using `n_mats=3`.
- Lines 9: computes `max_dim` using `max_dim=32`.
- Lines 10: computes `n_stats` using `n_stats=100`.
- Lines 11: computes `nnz_per_col` using `nnz_per_col=5`.
- Lines 16: computes `runtimes_full_poly` using `runtimes_full_poly=zeros(n_stats,1)`.
- Lines 17: computes `runtimes_full_flat` using `runtimes_full_flat=zeros(n_stats,1)`.
- Lines 26: computes `terms` using `terms=cell(1,n_mats)`.
- Lines 28: computes `dim` using `dim=randi(max_dim,1)`.
- Lines 29-30: computes `terms{k}` using `terms{k}=1i*randn(dim,dim)+ randn(dim,dim)`.
- Lines 34: computes `P` using `P=polyadic({terms})`.
- Lines 37-38: computes `v` using `v=1i*randn(size(P,2),1)+ randn(size(P,2),1)`.
- Lines 41: computes `tic; P*v; runtimes_full_poly(n)` using `tic; P*v; runtimes_full_poly(n)=1000*toc`.
- Lines 47: computes `tic; P*v; runtimes_full_flat(n)` using `tic; P*v; runtimes_full_flat(n)=1000*toc`.
- Lines 54: computes `runtimes_sparse_poly` using `runtimes_sparse_poly=zeros(n_stats,1)`.
- Lines 55: computes `runtimes_sparse_flat` using `runtimes_sparse_flat=zeros(n_stats,1)`.
- Lines 78: computes `tic; P*v; runtimes_sparse_poly(n)` using `tic; P*v; runtimes_sparse_poly(n)=1000*toc`.
- Lines 84: computes `tic; P*v; runtimes_sparse_flat(n)` using `tic; P*v; runtimes_sparse_flat(n)=1000*toc`.

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
