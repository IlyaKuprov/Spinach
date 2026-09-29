# examples/benchmarks/parallelization_1.m

## Use

`parallelization_1()` takes no arguments and returns no values. It is a timing example for evaluating observables during Hilbert-space time propagation on different local parallel-pool sizes. It requires the Spinach setup used by the example and MATLAB parallel-pool support for the requested sizes.

## Spin system and propagation

The example defines 12 1H spins for 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one, with sys.magnet=14.095, 12 assigned inter.zeeman.scalar values ([8.345, 7.741, 8.097, 8.354, 7.784, 8.330, 7.059, 7.941, 7.466, 7.326, 7.466, 7.941]), and the following explicitly assigned scalar-coupling entries:

| Spin pair | Value |
|---|---:|
| (1,2), (2,3), (8,9), (9,10), (10,11), (11,12) | 7.8 |
| (1,3) | 0.9 |
| (4,5) | 8.4 |
| (4,6), (8,10), (10,12) | 1.2 |
| (5,6) | 7.2 |

The source does not annotate units for these numeric settings, so the values above are reproduced without unit labels. It selects zeeman-hilb with none approximation, applies the nmr assumption, builds the Hamiltonian, and uses operator(spin_system,'Lx','all') as both the initial state and observable.

## Parallel timing

The candidate pool sizes are 2, 4, 8, 16, 32, 64, 128, 256, 512, and 1,024, filtered to values no greater than feature('numcores'). For each retained size it deletes the current pool, starts a pool of that size, waits 10 seconds, then times one call to evolution(spin_system,H,rho,rho,1e-3,1000,'observable') and prints the elapsed seconds. Pool creation and the explicit 10-second pause are outside the timed interval. The example replaces any already-running local pool and leaves the last created pool open when the loop completes. It filters by the reported core count, not by a pool profile's resource or licensing limits; parpool(n) must still be available for the selected size.

## Reference and source

The source cites Penchav et al., Spectrochimica Acta Part A 78 (2011), 559–565, [doi:10.1063/1.3679656](https://doi.org/10.1063/1.3679656).

[examples/benchmarks/parallelization_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/parallelization_1.m)
