# examples/fundamentals/quadratures/mas_benchmark.m

- Signature: `mas_benchmark()`

## Purpose

Compares propagation schemes over one period of a 50 kHz magic-angle-spinning rotor, using a fine-grid propagation as the reference.

## Physical / mathematical content

- The example models two `1H` spins at 14.1 T, with Zeeman shifts 5.0 and -2.0 and coordinates `[0 0 0]` and `[0 3.9 0.1]`. It builds the `sphten-liouv` basis with no approximation and projection +1, then applies the NMR assumption.
- The rotor period is `T=1/50000` s and the magic angle is `atan(sqrt(2))`. The initial state is `L+` on `1H`.

## Numerical / algorithmic content

- Hamiltonians are sampled at rotor phases across the period. A reference trajectory uses `2^13+1` samples and three successive Hamiltonians per `step` call.
- For grids `2^n+1`, `n=4:12`, the code measures relative state error for the left-point and midpoint piecewise-constant schemes and second- and fourth-order Lie-group schemes (LP, MP, LG-2, LG-4).

## Implementation structure

- Constructs the spin system and Hamiltonian, precomputes the rotor-oriented Hamiltonian samples, propagates the reference and benchmark trajectories, then plots relative error against the number of time-grid points on log-log axes.
- The source estimates a runtime of seconds.
