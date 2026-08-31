# examples/nmr_paramag/carb_anh/s217c_kuprov.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/carb_anh/s217c_kuprov.m`
- Signature: `s217c_kuprov()`
- Total lines: 52

## Purpose

Distributed fit for the S217C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Load experimental data; implemented by `load('s217c_expt.mat','expt_pcs','xyz','xyz_all')`.
- Lines 23-24: Load susceptibility tensor; implemented by `load('s217c_chi_eff.mat','chi')`.
- Lines 26-27: Set inverse problem parameters; implemented by `parameters.equation='kuprov'`.
- Lines 41-42: Solve and refine the grid; implemented by `for n=[64 128 256]`.
- Lines 47-48: Get the new susceptibility tensor; implemented by `[chi,~]=chi_eff(source_cube,ranges,xyz,expt_pcs)`.

### Control flow inferred from the code

- Line 42: `for` loop over `n=[64 128 256]`.

### Key state/data transformations

- Lines 27: computes `parameters.equation` using `parameters.equation='kuprov'`.
- Lines 28-29: computes `parameters.plot` using `parameters.plot={'diagnostics','density', 'molecule','tightzoom','box'}`.
- Lines 30: computes `parameters.box_cent` using `parameters.box_cent=[-21.8 -18.4 20.2]`.
- Lines 31: computes `parameters.box_size` using `parameters.box_size=[ 25.0 25.0 25.0]`.
- Lines 32: computes `parameters.margins` using `parameters.margins=50*ones(1,6)`.
- Lines 33: computes `parameters.confine` using `parameters.confine=[3.0 12.0]`.
- Lines 34: computes `parameters.sharpen` using `parameters.sharpen=1.0`.
- Lines 35: computes `parameters.xyz_all` using `parameters.xyz_all=xyz_all`.
- Lines 36: computes `parameters.expt_pcs` using `parameters.expt_pcs=expt_pcs`.
- Lines 37: computes `parameters.xyz` using `parameters.xyz=xyz`.
- Lines 38: computes `parameters.chi` using `parameters.chi=chi`.
- Lines 39: computes `parameters.gpu` using `parameters.gpu=true()`.
- Lines 43: computes `[source_cube,ranges]` using `[source_cube,ranges]=ipcs(parameters,n,0.50)`.
- Lines 44: computes `parameters.guess` using `parameters.guess=source_cube`.
- Lines 48: computes `[chi,~]` using `[chi,~]=chi_eff(source_cube,ranges,xyz,expt_pcs)`.

## Implementation structure

- Distributed fit for the S217C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Load susceptibility tensor
- Set inverse problem parameters
- Solve and refine the grid
- Get the new susceptibility tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `true()`, `ipcs()`, `chi_eff()`.
