# examples/esr_sol_swept/fieldsweep_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_nitroxide.m`
- Signature: `fieldsweep_nitroxide()`
- Total lines: 56

## Purpose

Field swept EPR spectrum of nitroxide, computed by finding resonance fields and transition moments. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Field swept EPR spectrum of nitroxide, computed by finding
- resonance fields and transition moments.
- Calculation time: seconds.
- Isotopes
- Magnet field (must be 1)
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the simulation in the high-T approximation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
