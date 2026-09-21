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
