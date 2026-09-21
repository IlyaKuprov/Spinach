# examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/calbindin/tm_1igv_mult_fit.m`
- Signature: `tm_1igv_mult_fit()`
- Total lines: 29

## Purpose

Electron location and susceptibility tensor recovery from experimental PCS data using point electron model. Experi- mental data kindly provided by Gottfried Otting (Australi- an National University).

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

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
