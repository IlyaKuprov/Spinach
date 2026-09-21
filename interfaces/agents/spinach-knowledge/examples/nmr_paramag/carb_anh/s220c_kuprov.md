# examples/nmr_paramag/carb_anh/s220c_kuprov.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/carb_anh/s220c_kuprov.m`
- Signature: `s220c_kuprov()`
- Total lines: 52

## Purpose

Distributed fit for the S220C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

## Implementation structure

- Distributed fit for the S220C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Load susceptibility tensor
- Set inverse problem parameters
- Solve and refine the grid
- Get the new susceptibility tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `true()`, `ipcs()`, `chi_eff()`.
