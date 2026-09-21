# examples/esr_sol_swept/fieldsweep_triplet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_triplet.m`
- Signature: `fieldsweep_triplet()`
- Total lines: 64

## Purpose

Powder averaged X-band field-swept ESR spectrum of photo- generated pentacene triplet state. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder averaged X-band field-swept ESR spectrum of photo-
- generated pentacene triplet state.
- Calculation time: seconds.
- Magnet field (must be 1)
- Triplet electron
- Zeeman tensor, assumed isotropic
- ZFS, photo-excited pentacene triplet
- Basis set
- Spinach housekeeping
- Experiment parameters
- Zeeman tensor into Hz/Tesla
- Orientation-and field-dependent initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs2mat()`, `create()`, `basis()`, `spin()`, `zftrip()`, `euler2dcm()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
