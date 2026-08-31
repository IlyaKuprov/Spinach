# examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m`
- Signature: `vfit_4sim_ik()`
- Total lines: 168

## Purpose

Simultaneous fitting of multiple 51V MAS NMR spectra with respect to the chemical shielding anisotropy and quadrupole coupling tensor parameters. Calculation time: hours, much faster with a GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Load and filter the data; implemented by `s29=load('v12_29_dec15.spc'); s29=sgolayfilt(s29,3,51)`.
- Lines 18-19: Set spectral ranges; implemented by `kstart = 6200`.
- Lines 27-28: Preprocess the spectra; implemented by `S35=zeros(nrpoints_b,1); A35=s35(myrangeb,1)`.
- Lines 37-38: Set the initial guess; implemented by `guess=[-669.0 564.0 0.255 82.0 180.0 19.0 3.72 0.62]`.
- Lines 40-41: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 43-44: Get a figure going; implemented by `kfigure(); scale_figure([1.75 1.50])`.
- Lines 46-47: Run the optimisation; implemented by `fminsearch(@(x)errfun(A35,A33,A31,A29,S35,S33,S31,S29,x),guess,options)`.

### Key state/data transformations

- Lines 13: computes `s29` using `s29=load('v12_29_dec15.spc'); s29=sgolayfilt(s29,3,51)`.
- Lines 14: computes `s31` using `s31=load('v12_31_dec15.spc'); s31=sgolayfilt(s31,3,51)`.
- Lines 15: computes `s33` using `s33=load('v12_33_dec15.spc'); s33=sgolayfilt(s33,3,51)`.
- Lines 16: computes `s35` using `s35=load('v12_35_dec15.spc'); s35=sgolayfilt(s35,3,51)`.
- Lines 19: computes `kstart` using `kstart = 6200`.
- Lines 20: computes `kend_a` using `kend_a = 10000`.
- Lines 21: computes `kend_b` using `kend_b = kstart+4096-1`.
- Lines 22: computes `sp_range` using `sp_range = kstart:kend_a`.
- Lines 23: computes `myrangeb` using `myrangeb = kstart:kend_b`.
- Lines 24: computes `nrpoints_a` using `nrpoints_a = kend_a-kstart+1`.
- Lines 25: computes `nrpoints_b` using `nrpoints_b = 4096`.
- Lines 28: computes `S35` using `S35=zeros(nrpoints_b,1); A35=s35(myrangeb,1)`.
- Lines 29: computes `S33` using `S33=zeros(nrpoints_b,1); A33=s33(myrangeb,1)`.
- Lines 30: computes `S31` using `S31=zeros(nrpoints_b,1); A31=s31(myrangeb,1)`.
- Lines 31: computes `S29` using `S29=zeros(nrpoints_b,1); A29=s29(myrangeb,1)`.
- Lines 32: computes `S35(1:nrpoints_a,1)` using `S35(1:nrpoints_a,1)=s35(sp_range,2)/max(s35(sp_range,2))`.
- Lines 33: computes `S33(1:nrpoints_a,1)` using `S33(1:nrpoints_a,1)=s33(sp_range,2)/max(s33(sp_range,2))`.
- Lines 34: computes `S31(1:nrpoints_a,1)` using `S31(1:nrpoints_a,1)=s31(sp_range,2)/max(s31(sp_range,2))`.

### Local helper functions

- Line 53: `errfun()` — `function err=errfun(A35,A33,A31,A29,S35,S33,S31,S29,params)`. Display the parameters
  - Representative operation: `disp(params)`.
  - Representative operation: `sys.output='hush'`.

## Implementation structure

- Simultaneous fitting of multiple 51V MAS NMR spectra with
- respect to the chemical shielding anisotropy and quadrupole
- coupling tensor parameters.
- Calculation time: hours, much faster with a GPU.
- Load and filter the data
- Set spectral ranges
- Preprocess the spectra
- Set the initial guess
- Set optimiser options
- Get a figure going
- Run the optimisation
- Least squares error function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `sgolayfilt()`, `s35()`, `s33()`, `s31()`, `s29()`, `S35()`, `S33()`, `S31()`, `S29()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `A35()`.
