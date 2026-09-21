# examples/esr_sol_swept/fieldsweep_porphyrin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_porphyrin.m`
- Signature: `fieldsweep_porphyrin()`
- Total lines: 69

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
