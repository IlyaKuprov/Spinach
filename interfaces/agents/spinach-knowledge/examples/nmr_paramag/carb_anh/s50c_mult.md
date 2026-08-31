# examples/nmr_paramag/carb_anh/s50c_mult.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/carb_anh/s50c_mult.m`
- Signature: `s50c_mult()`
- Total lines: 38

## Purpose

Multipolar fit for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Load experimental data; implemented by `load('s50c_expt.mat','expt_pcs','xyz')`.
- Lines 23-24: Solve the inverse problem; implemented by `[mxyz,chi,~,pred_pcs]=ilpcs(xyz,expt_pcs,[0 1 2],[-27.0 13.0 18.0])`.
- Lines 26-27: Plot experimental vs predicted PCS; implemented by `kfigure(); plot(expt_pcs,pred_pcs,'bo'); hold on; kgrid`.
- Lines 33-34: Report and save the parameters; implemented by `disp('Susceptibility tensor:'); disp(chi)`.

### Key state/data transformations

- Lines 24: computes `[mxyz,chi,~,pred_pcs]` using `[mxyz,chi,~,pred_pcs]=ilpcs(xyz,expt_pcs,[0 1 2],[-27.0 13.0 18.0])`.

## Implementation structure

- Multipolar fit for the S50C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Solve the inverse problem
- Plot experimental vs predicted PCS
- Report and save the parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `ilpcs()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`, `ylim()`.
