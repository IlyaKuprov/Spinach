# examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m`
- Signature: `kbr_mas_fitting()`
- Total lines: 113

## Purpose

Fitting of a 79Br MAS NMR spectrum of potassium bromide with respect to the quadrupole coupling constant. The spectrum cannot be fitted with a single quadrupolar tensor; at least 3 are necessary, likely due to a dist- ribution of electrostatic environments in the powder. Calculation time: hours.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Load and normalise the data; implemented by `kbr_data=readmatrix('KBr_400MHz_2kHz.txt')`.
- Lines 20-21: Set instrumental variables; implemented by `sys.magnet=9.3659`.
- Lines 35-36: Set optimizer options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 38-39: Get a figure going; implemented by `kfigure(); scale_figure([1.75 1.50])`.
- Lines 41-42: Run the optimisation; implemented by `fminsearch(@errfun,guess,options)`.

### Key state/data transformations

- Lines 17: computes `kbr_data` using `kbr_data=readmatrix('KBr_400MHz_2kHz.txt')`.
- Lines 18: computes `spec_expt` using `spec_expt=kbr_data(:,2)/100`.
- Lines 21: computes `sys.magnet` using `sys.magnet=9.3659`.
- Lines 22: computes `parameters.rate` using `parameters.rate=2000`.
- Lines 23: computes `parameters.offset` using `parameters.offset=6034.96`.
- Lines 24: computes `parameters.sweep` using `parameters.sweep=1e5`.
- Lines 25: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 26: computes `parameters.zerofill` using `parameters.zerofill=32768`.
- Lines 28: computes `guess` using `guess=[60.0933`.
- Lines 36: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.

### Local helper functions

- Line 45: `errfun()` — `function err=errfun(params)`. Silence Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene','trajlevel'}`.

## Implementation structure

- Fitting of a 79Br MAS NMR spectrum of potassium bromide
- with respect to the quadrupole coupling constant.
- The spectrum cannot be fitted with a single quadrupolar
- tensor; at least 3 are necessary, likely due to a dist-
- ribution of electrostatic environments in the powder.
- Calculation time: hours.
- Load and normalise the data
- Set instrumental variables
- Set optimizer options
- Get a figure going
- Run the optimisation
- Least squares error function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `readmatrix()`, `kbr_data()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `plot_1d()`.
