# examples/kinetics/frydman_pump_b.m

- Signature: `frydman_pump_b()`

## Purpose

Lucio Frydman's water-exchange spin-lock pump, Figure 9 from https://doi.org/10.1016/j.jmr.2021.107083. The script simulates ten repeated pump cycles and plots peptide H, N, and carbonyl-C Z/X trajectories; the source estimates seconds.

## Physical / mathematical content

The spin system contains peptide H, N, Cα, and C′ plus 100 water protons. The source assigns H–N and N–C′ couplings of −45 and 8 Hz, respectively, and builds intermolecular exchange with NH–water and water-pool rates of 1000 and 10⁴ Hz. Relaxation uses the specified diagonal T1/T2 rates and isotropic equilibrium at 298 K.

## Numerical / algorithmic content

Each cycle applies 90° forward flips, two CP periods (11 ms and 53 ms), reverse flips, and a crusher. The two CP stages use 11 and 53 points. Separate effective Hamiltonians are formed for the CP periods, and `evolution` generates their trajectories.

## Implementation structure

- Uses B₀ = 11.7 T, the H/N/C channels, and a 100-proton water pool.
- Sets `parameters.nloops=10`; CP durations are [11, 53] ms with [11, 53] points.
- Computes and plots H, N, and carbonyl-C Z/X expectation values against time.
