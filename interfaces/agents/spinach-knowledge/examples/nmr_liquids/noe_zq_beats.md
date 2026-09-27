# examples/nmr_liquids/noe_zq_beats.m

- Signature: `noe_zq_beats()`

## Purpose

Zero-quantum beats in the Overhauser effect in a strongly coupled two-spin system. Calculation time: seconds.

## Physical / mathematical content

- The system is two protons with a 0.01 ppm difference in isotropic shift, a 3.0 Hz scalar coupling, and a 2.00 Å separation. Redfield relaxation is specified with a 1 ns correlation time, 298 K, the Di Bari equilibrium convention, and secular retention.
- Starting from thermal equilibrium with spin 1 inverted, the calculation follows both longitudinal magnetizations; the strong coupling allows the zero-quantum-beat behaviour associated with the NOE to appear.

## Numerical / algorithmic content

- The Liouvillian is the NMR-frame Hamiltonian plus the relaxation contribution. It is propagated in the full spherical-tensor Liouville basis, with a 4.0 proximity cutoff.
- The two longitudinal detection channels are sampled every 0.01 s over 1000 intervals (0–10 s).

## Implementation structure

- Build the two-proton system at 14.1 T, construct the Redfield superoperator, and add it to the assumed NMR Hamiltonian.
- Invert spin 1 relative to equilibrium, run multichannel evolution with both `Lz` operators, and plot the real longitudinal signals for Proton A and Proton B.
