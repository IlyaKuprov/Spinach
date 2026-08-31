# examples/nmr_paramag/carb_anh/s50c_lcurve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/carb_anh/s50c_lcurve.m`
- Signature: `s50c_lcurve()`
- Total lines: 57

## Purpose

L-curves for the S50C mutant dataset for human carbonic anhydrase II. The system and the method are described in: A step-by-step tutorial is available here:

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Load experimental data; implemented by `load('s50c_expt.mat','expt_pcs','xyz','xyz_all')`.
- Lines 23-24: Load susceptibility tensor; implemented by `load('s50c_chi_eff.mat','chi')`.
- Lines 26-27: Solver parameters; implemented by `parameters.plot={}`.
- Lines 40-41: Regularisation parameter array; implemented by `lam=10.^linspace(-2.0,2.0,30)`.
- Lines 43-44: Result arrays; implemented by `err=zeros(1,30); reg=zeros(1,30)`.
- Lines 46-47: Run a parallel loop; implemented by `parfor n=1:30`.
- Lines 52-53: L-curve analysis; implemented by `kfigure(); S=lcurve(lam,err,reg,'log'); drawnow`.

### Control flow inferred from the code

- Line 47: `parfor` loop over `n=1:30`.

### Key state/data transformations

- Lines 27: computes `parameters.plot` using `parameters.plot={}`.
- Lines 28: computes `parameters.equation` using `parameters.equation='kuprov'`.
- Lines 29: computes `parameters.box_cent` using `parameters.box_cent=[-27.4 13.3 18.8]`.
- Lines 30: computes `parameters.box_size` using `parameters.box_size=[ 50.0 50.0 50.0]`.
- Lines 31: computes `parameters.margins` using `parameters.margins=50*ones(1,6)`.
- Lines 32: computes `parameters.confine` using `parameters.confine=[2.0 12.0]`.
- Lines 33: computes `parameters.sharpen` using `parameters.sharpen=0.0`.
- Lines 34: computes `parameters.xyz_all` using `parameters.xyz_all=xyz_all`.
- Lines 35: computes `parameters.expt_pcs` using `parameters.expt_pcs=expt_pcs`.
- Lines 36: computes `parameters.xyz` using `parameters.xyz=xyz`.
- Lines 37: computes `parameters.chi` using `parameters.chi=chi`.
- Lines 38: computes `parameters.gpu` using `parameters.gpu=true()`.
- Lines 41: computes `lam` using `lam=10.^linspace(-2.0,2.0,30)`.
- Lines 44: computes `err` using `err=zeros(1,30); reg=zeros(1,30)`.
- Lines 48: computes `[~,~,~,err(n),~,reg(n)]` using `[~,~,~,err(n),~,reg(n)]=ipcs(parameters,64,lam(n))`.
- Lines 49: computes `reg(n)` using `reg(n)=reg(n)/lam(n)`.
- Lines 53: computes `kfigure(); S` using `kfigure(); S=lcurve(lam,err,reg,'log'); drawnow`.

## Implementation structure

- L-curves for the S50C mutant dataset for human carbonic anhydrase
- II. The system and the method are described in:
- A step-by-step tutorial is available here:
- Load experimental data
- Load susceptibility tensor
- Solver parameters
- Regularisation parameter array
- Result arrays
- Run a parallel loop
- L-curve analysis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `true()`, `err()`, `reg()`, `ipcs()`, `lam()`, `kfigure()`, `lcurve()`, `num2str()`.
