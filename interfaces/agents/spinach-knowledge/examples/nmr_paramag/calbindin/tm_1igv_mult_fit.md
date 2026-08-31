# examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m`
- Signature: `tm_1igv_mult_fit()`
- Total lines: 29

## Purpose

Electron location and susceptibility tensor recovery from experimental PCS data using point electron model. Experi- mental data kindly provided by Gottfried Otting (Australi- an National University).

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Load experimental data; implemented by `load('tm_1igv_pcs.mat','expt_pcs','x','y','z')`.
- Lines 14-15: Solve the inverse problem; implemented by `[mxyz,chi,~,pred_pcs]=ilpcs([x y z],expt_pcs,[0 1 2],[-5 5 -15])`.
- Lines 17-18: Plot experimental vs predicted PCS; implemented by `kfigure(); plot(expt_pcs,pred_pcs,'bo'); hold on; kgrid`.
- Lines 24-25: Report the parameters; implemented by `disp('Susceptibility tensor:'); disp(chi)`.

### Key state/data transformations

- Lines 15: computes `[mxyz,chi,~,pred_pcs]` using `[mxyz,chi,~,pred_pcs]=ilpcs([x y z],expt_pcs,[0 1 2],[-5 5 -15])`.

## Implementation structure

- Electron location and susceptibility tensor recovery from
- experimental PCS data using point electron model. Experi-
- mental data kindly provided by Gottfried Otting (Australi-
- an National University).
- Load experimental data
- Solve the inverse problem
- Plot experimental vs predicted PCS
- Report the parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `ilpcs()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`, `ylim()`.
