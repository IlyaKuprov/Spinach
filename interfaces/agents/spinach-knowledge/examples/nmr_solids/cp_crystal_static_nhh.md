# examples/nmr_solids/cp_crystal_static_nhh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_crystal_static_nhh.m`
- Signature: `cp_crystal_static_nhh()`
- Total lines: 67

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Static single crystal simulation in a full Liouvil- le space (here necessary because this is not a powder and everything interacts with everything). Calculation time: minutes on a Tesla A100, much longer on CPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
- sphere around it. Static single crystal simulation in a full Liouvil-
- le space (here necessary because this is not a powder and everything
- interacts with everything).
- Calculation time: minutes on a Tesla A100, much longer on CPU.
- System specification
- Interactions
- Basis set
- This needs a GPU
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `crystal()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
