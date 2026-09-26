# examples/esr_sol_pulsed/oop_eseem_photosystem_1.m

- Signature: `oop_eseem_photosystem_1()`

## Purpose

Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-] spin-correlated electron pair in Photosystem I. Time-domain si- mulation in Liouville space with averaging over a powder grid. Set to reproduce Figure 3a in The eigenvalues of P700+ g-tensor (without angles) come from The eigenvales of A-g-tensor (without angles) come from Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder-averaged two-pulse out-of-phase ESEEM on the [P700+,A1-]
- spin-correlated electron pair in Photosystem I. Time-domain si-
- mulation in Liouville space with averaging over a powder grid.
- Set to reproduce Figure 3a in
- The eigenvalues of P700+ g-tensor (without angles) come from
- The eigenvales of A-g-tensor (without angles) come from
- Calculation time: seconds
- Magnet field
- System specification
- Relaxation theory
- Basis set
- Disable trajectory-level SSR algorithms
