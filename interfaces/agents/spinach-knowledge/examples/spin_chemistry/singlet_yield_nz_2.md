# examples/spin_chemistry/singlet_yield_nz_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_nz_2.m`
- Signature: `singlet_yield_nz_2()`
- Total lines: 112

## Purpose

Field dependence of the decay rate of a micelle-confined triplet-born benzophenone ketyl / alkyl radical pair, computed with Redfield theory and with the lifetime-shifted Nakajima-Zwanzig kernel. The high-field decay of micellar pairs is relaxation-controlled: the T+/-states drain into the reactive S/T0 subspace at rates set by spectral densities of the anisotropic hyperfine modulation. Contact recombination at 1e9 H

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file also defines local helper function(s): `scrp_decay()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Field grid, contact recombination rate, and correlation time; implemented by `fields=[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.7 1.0 1.34]`.
- Lines 27-28: Preallocate observed decay rates; implemented by `decay_rf=zeros(size(fields))`.
- Lines 31-32: Loop over the field grid; implemented by `for n=1:numel(fields)`.
- Lines 34-35: Compute both theories at this field; implemented by `decay_rf(n)=scrp_decay(fields(n),k_rec,tau_c,'redfield')`.
- Lines 40-41: Report the two curves; implemented by `disp('field (T), Redfield decay rate (Hz), Nakajima-Zwanzig decay rate (Hz):')`.
- Lines 44-45: Plot the two curves; implemented by `kfigure(); semilogy(fields,[decay_rf' decay_nz'],'-o')`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=1:numel(fields)`.

### Key state/data transformations

- Lines 24: computes `fields` using `fields=[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.7 1.0 1.34]`.
- Lines 25: computes `k_rec` using `k_rec=1e9; tau_c=0.7e-9`.
- Lines 28: computes `decay_rf` using `decay_rf=zeros(size(fields))`.
- Lines 29: computes `decay_nz` using `decay_nz=zeros(size(fields))`.
- Lines 35: computes `decay_rf(n)` using `decay_rf(n)=scrp_decay(fields(n),k_rec,tau_c,'redfield')`.
- Lines 36: computes `decay_nz(n)` using `decay_nz(n)=scrp_decay(fields(n),k_rec,tau_c,'naka-zwan')`.

### Local helper functions

- Line 53: `scrp_decay()` — `function decay=scrp_decay(field,k_rec,tau_c,theory)`. Magnet and isotopes
  - Representative operation: `sys.magnet=field`.
  - Representative operation: `sys.isotopes={'E','E','1H','1H'}`.

## Implementation structure

- Field dependence of the decay rate of a micelle-confined triplet-born
- benzophenone ketyl / alkyl radical pair, computed with Redfield theory
- and with the lifetime-shifted Nakajima-Zwanzig kernel. The high-field
- decay of micellar pairs is relaxation-controlled: the T+/-states drain
- into the reactive S/T0 subspace at rates set by spectral densities of
- the anisotropic hyperfine modulation. Contact recombination at 1e9 Hz
- against a supercage correlation time of 0.7 ns puts the scalar lifetime
- shift at k*tau_c near 0.35, and the two theories separate. The observed
- decay rate is the slowest eigenmode of the pair Liouvillian. Parameters
- are representative of the SDS supercage systems of Sakaguchi, Hayashi,
- and Nagakura:
- Calculation time: minutes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `decay_rf()`, `scrp_decay()`, `fields()`, `decay_nz()`, `field()`, `rate()`, `num2str()`, `kfigure()`, `semilogy()`, `klegend()`, `kxlabel()`, `kylabel()`, `mt2hz()`, `strcmp()`, `create()`, `basis()`.
