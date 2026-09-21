# examples/nmr_liquids/noe_two_spin_hom.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noe_two_spin_hom.m`
- Signature: `noe_two_spin_hom()`
- Total lines: 56

## Purpose

Nuclear overhauser effect in a homonuclear two-spin system in the long correlation time case. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Nuclear overhauser effect in a homonuclear two-spin system in
- the long correlation time case.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Spinach housekeeping
- Build the relaxation superoperator
- Get thermal equilibrium state
- Start in a state with one spin inverted
- Compute the evolution trajectory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`, `evolution()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
