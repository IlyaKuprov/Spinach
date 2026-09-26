# examples/esr_sol_swept/fieldsweep_porphyrin.m

- Signature: `fieldsweep_porphyrin()`

## Purpose

Field swept EPR spectrum of copper porphyrin complex, computed by finding resonance fields and transition moments. Calculation time: minutes.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Field swept EPR spectrum of copper porphyrin complex, computed
- by finding resonance fields and transition moments.
- Calculation time: minutes.
- Magnet field
- Isotopes
- Array preallocation
- Zeeman interactions
- Hyperfine interactions
- Basis set
- Symmetry
- Spinach housekeeping
- Experiment parameters
