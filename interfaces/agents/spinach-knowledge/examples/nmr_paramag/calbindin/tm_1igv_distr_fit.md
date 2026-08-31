# examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/calbindin/tm_1igv_distr_fit.m`
- Signature: `tm_1igv_distr_fit()`
- Total lines: 48

## Purpose

Inverse problem for the unpaired electron density distribution. Experimental data from Gottfried Ott- ing (Australian National University).

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Load the pdb file; implemented by `pdb=pdbread('1igv_processed.pdb')`.
- Lines 13-14: Load experimental data; implemented by `load('tm_1igv_pcs.mat','x','y','z','expt_pcs')`.
- Lines 16-17: Load susceptibility tensor; implemented by `load('tm_1igv_chi_eff.mat','chi')`.
- Lines 19-20: Inverse solver parameters; implemented by `parameters.equation='kuprov'`.
- Lines 36-37: Iteratively refine the grid; implemented by `for n=[64 128 256 384]`.
- Lines 42-43: Get the new susceptibility tensor; implemented by `[chi,~]=chi_eff(source_cube,ranges,[x y z],expt_pcs)`.

### Control flow inferred from the code

- Line 37: `for` loop over `n=[64 128 256 384]`.

### Key state/data transformations

- Lines 11: computes `pdb` using `pdb=pdbread('1igv_processed.pdb')`.
- Lines 20: computes `parameters.equation` using `parameters.equation='kuprov'`.
- Lines 21-22: computes `parameters.plot` using `parameters.plot={'diagnostics','density', 'molecule','tightzoom','box'}`.
- Lines 23: computes `parameters.box_cent` using `parameters.box_cent=[3.5 17.0 16.1]`.
- Lines 24: computes `parameters.box_size` using `parameters.box_size=[7.0 7.0 7.0]`.
- Lines 25-27: computes `parameters.xyz_all` using `parameters.xyz_all=[[pdb.Model.Atom(:).X]' [pdb.Model.Atom(:).Y]' [pdb.Model.Atom(:).Z]']`.
- Lines 28: computes `parameters.margins` using `parameters.margins=50*ones(1,6)`.
- Lines 29: computes `parameters.confine` using `parameters.confine=[1.0 3.0]`.
- Lines 30: computes `parameters.sharpen` using `parameters.sharpen=2e3`.
- Lines 31: computes `parameters.expt_pcs` using `parameters.expt_pcs=expt_pcs`.
- Lines 32: computes `parameters.xyz` using `parameters.xyz=[x y z]`.
- Lines 33: computes `parameters.chi` using `parameters.chi=chi`.
- Lines 34: computes `parameters.gpu` using `parameters.gpu=true()`.
- Lines 38: computes `[source_cube,ranges]` using `[source_cube,ranges]=ipcs(parameters,n,0.34)`.
- Lines 39: computes `parameters.guess` using `parameters.guess=source_cube`.
- Lines 43: computes `[chi,~]` using `[chi,~]=chi_eff(source_cube,ranges,[x y z],expt_pcs)`.

## Implementation structure

- Inverse problem for the unpaired electron density
- distribution. Experimental data from Gottfried Ott-
- ing (Australian National University).
- Load the pdb file
- Load experimental data
- Load susceptibility tensor
- Inverse solver parameters
- Iteratively refine the grid
- Get the new susceptibility tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pdbread()`, `load()`, `true()`, `ipcs()`, `chi_eff()`, `save()`.
