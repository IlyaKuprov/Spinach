# examples/nmr_paramag/calbindin/tm_1igv_lcurve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/calbindin/tm_1igv_lcurve.m`
- Signature: `tm_1igv_lcurve()`
- Total lines: 52

## Purpose

Inverse problem for the unpaired electron density distribution. Experimental data from Gottfried Ott- ing (Australian National University).

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Load the pdb file; implemented by `pdb=pdbread('1igv_processed.pdb')`.
- Lines 13-14: Load experimental data; implemented by `load('tm_1igv_pcs.mat','x','y','z','expt_pcs')`.
- Lines 16-17: Load susceptibility tensor; implemented by `load('tm_1igv_chi_eff.mat','chi')`.
- Lines 19-20: Inverse solver parameters; implemented by `parameters.plot={}`.
- Lines 35-36: Regularisation parameter array; implemented by `lam=10.^linspace(-1.5,1.0,15)`.
- Lines 38-39: Result arrays; implemented by `err=zeros(1,15); reg=zeros(1,15)`.
- Lines 41-42: Run a parallel loop; implemented by `parfor n=1:15`.
- Lines 47-48: L-curve analysis; implemented by `S=lcurve(lam,err,reg,'log'); drawnow`.

### Control flow inferred from the code

- Line 42: `parfor` loop over `n=1:15`.

### Key state/data transformations

- Lines 11: computes `pdb` using `pdb=pdbread('1igv_processed.pdb')`.
- Lines 20: computes `parameters.plot` using `parameters.plot={}`.
- Lines 21: computes `parameters.equation` using `parameters.equation='kuprov'`.
- Lines 22: computes `parameters.box_cent` using `parameters.box_cent=[3.5 17.0 16.1]`.
- Lines 23: computes `parameters.box_size` using `parameters.box_size=[7.0 7.0 7.0]`.
- Lines 24-26: computes `parameters.xyz_all` using `parameters.xyz_all=[[pdb.Model.Atom(:).X]' [pdb.Model.Atom(:).Y]' [pdb.Model.Atom(:).Z]']`.
- Lines 27: computes `parameters.margins` using `parameters.margins=50*ones(1,6)`.
- Lines 28: computes `parameters.confine` using `parameters.confine=[1.0 3.0]`.
- Lines 29: computes `parameters.sharpen` using `parameters.sharpen=0.0`.
- Lines 30: computes `parameters.expt_pcs` using `parameters.expt_pcs=expt_pcs`.
- Lines 31: computes `parameters.xyz` using `parameters.xyz=[x y z]`.
- Lines 32: computes `parameters.chi` using `parameters.chi=chi`.
- Lines 33: computes `parameters.gpu` using `parameters.gpu=true()`.
- Lines 36: computes `lam` using `lam=10.^linspace(-1.5,1.0,15)`.
- Lines 39: computes `err` using `err=zeros(1,15); reg=zeros(1,15)`.
- Lines 43: computes `[~,~,~,err(n),~,reg(n)]` using `[~,~,~,err(n),~,reg(n)]=ipcs(parameters,128,lam(n))`.
- Lines 44: computes `reg(n)` using `reg(n)=reg(n)/lam(n)`.
- Lines 48: computes `S` using `S=lcurve(lam,err,reg,'log'); drawnow`.

## Implementation structure

- Inverse problem for the unpaired electron density
- distribution. Experimental data from Gottfried Ott-
- ing (Australian National University).
- Load the pdb file
- Load experimental data
- Load susceptibility tensor
- Inverse solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pdbread()`, `load()`, `true()`, `err()`, `reg()`, `ipcs()`, `lam()`, `lcurve()`, `num2str()`.
