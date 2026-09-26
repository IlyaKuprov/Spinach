# examples/nmr_liquids/noe_zq_beats.m

- Signature: `noe_zq_beats()`

## Purpose

Zero-quantum beats in the Overhauser effect in a strongly coupled two-spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Zero-quantum beats in the Overhauser effect in a strongly
- coupled two-spin system.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Build the Liouvillian
- Get thermal equilibrium state
- Start in a state with one spin inverted
