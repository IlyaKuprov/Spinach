# examples/spin_chemistry/singlet_yield_nz_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_nz_1.m`
- Signature: `singlet_yield_nz_1()`
- Total lines: 130

## Purpose

Magnetic field effect on a triplet-born benzophenone ketyl / thiyl radical pair in a viscous ionic liquid, computed with Redfield theory and with the lifetime-shifted Nakajima-Zwanzig kernel. The cage life- time of tens of nanoseconds and the rotational correlation time of a few nanoseconds put k*tau_c near 0.1, where the recombination drain visibly competes with rotational decorrelation of the anisotropic hyperfine 

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `cage_pair()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Field grid, cage drain, and correlation time; implemented by `fields=linspace(0.5e-3,50e-3,15)`.
- Lines 24-25: Preallocate singlet yields; implemented by `yield_rf=zeros(size(fields))`.
- Lines 28-29: Loop over the field grid; implemented by `for n=1:numel(fields)`.
- Lines 31-32: Compute both theories at this field; implemented by `yield_rf(n)=cage_pair(fields(n),k_cage,tau_c,'redfield')`.
- Lines 37-38: Report the two curves; implemented by `disp('field (T), Redfield yield, Nakajima-Zwanzig yield:')`.
- Lines 41-42: Drain rate ladder at the rising edge of the MARY curve; implemented by `drains=[1e6 3e6 1e7 3e7 1e8]`.
- Lines 49-50: Report the drain ladder; implemented by `disp('drain rate (Hz), Redfield yield, Nakajima-Zwanzig yield:')`.
- Lines 53-54: Plot the MARY curves and the drain ladder; implemented by `kfigure(); scale_figure([1.75 0.75])`.

### Control flow inferred from the code

- Line 29: `for` loop over `n=1:numel(fields)`.
- Line 44: `for` loop over `n=1:numel(drains)`.

### Key state/data transformations

- Lines 21: computes `fields` using `fields=linspace(0.5e-3,50e-3,15)`.
- Lines 22: computes `k_cage` using `k_cage=3e7; tau_c=3.3e-9`.
- Lines 25: computes `yield_rf` using `yield_rf=zeros(size(fields))`.
- Lines 26: computes `yield_nz` using `yield_nz=zeros(size(fields))`.
- Lines 32: computes `yield_rf(n)` using `yield_rf(n)=cage_pair(fields(n),k_cage,tau_c,'redfield')`.
- Lines 33: computes `yield_nz(n)` using `yield_nz(n)=cage_pair(fields(n),k_cage,tau_c,'naka-zwan')`.
- Lines 42: computes `drains` using `drains=[1e6 3e6 1e7 3e7 1e8]`.
- Lines 43: computes `edge_rf` using `edge_rf=zeros(size(drains)); edge_nz=zeros(size(drains))`.
- Lines 45: computes `edge_rf(n)` using `edge_rf(n)=cage_pair(3e-3,drains(n),tau_c,'redfield')`.
- Lines 46: computes `edge_nz(n)` using `edge_nz(n)=cage_pair(3e-3,drains(n),tau_c,'naka-zwan')`.

### Local helper functions

- Line 66: `cage_pair()` — `function yield=cage_pair(field,k_cage,tau_c,theory)`. Magnet and isotopes
  - Representative operation: `sys.magnet=field`.
  - Representative operation: `sys.isotopes={'E','E','1H','1H'}`.

## Implementation structure

- Magnetic field effect on a triplet-born benzophenone ketyl / thiyl
- radical pair in a viscous ionic liquid, computed with Redfield theory
- and with the lifetime-shifted Nakajima-Zwanzig kernel. The cage life-
- time of tens of nanoseconds and the rotational correlation time of a
- few nanoseconds put k*tau_c near 0.1, where the recombination drain
- visibly competes with rotational decorrelation of the anisotropic
- hyperfine couplings. Parameters are representative of the TMPA-TFSA
- measurements of the Wakasa group:
- Calculation time: minutes
- Field grid, cage drain, and correlation time
- Preallocate singlet yields
- Loop over the field grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `yield_rf()`, `cage_pair()`, `fields()`, `yield_nz()`, `field()`, `num2str()`, `edge_rf()`, `drains()`, `edge_nz()`, `rate()`, `kfigure()`, `scale_figure()`, `subplot()`, `klegend()`, `kxlabel()`, `kylabel()`.
