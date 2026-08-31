# examples/spin_chemistry/cidnp_nz_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_nz_1.m`
- Signature: `cidnp_nz_1()`
- Total lines: 124

## Purpose

Field dependence of geminate CIDNP from a radical pair in a viscous solvent, computed with Redfield theory and with the lifetime-shifted Nakajima-Zwanzig kernel. At a rotational correlation time of 1 ns and a singlet recombination rate of 1e8 Hz, the anisotropic hyperfine relaxation proceeds at a fair fraction of the recombination rate, and the pair drains before the bath decorrelates: the spectral densities seen by 

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `geminate_pol()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Field grid, recombination rate, and correlation time; implemented by `fields=[0.002 0.005 0.01 0.02 0.035 0.05 0.075]`.
- Lines 26-27: Preallocate product polarisations; implemented by `pol_rf=zeros(size(fields))`.
- Lines 30-31: Loop over the field grid; implemented by `for n=1:numel(fields)`.
- Lines 33-34: Compute both theories at this field; implemented by `pol_rf(n)=geminate_pol(fields(n),k_rec,tau_c,'redfield')`.
- Lines 39-40: Report the two curves; implemented by `disp('field (T), Redfield product polarisation, Nakajima-Zwanzig product polarisation:')`.
- Lines 43-44: Plot the two curves; implemented by `kfigure(); plot(1e3*fields,[pol_rf' pol_nz'],'-o')`.

### Control flow inferred from the code

- Line 31: `for` loop over `n=1:numel(fields)`.

### Key state/data transformations

- Lines 23: computes `fields` using `fields=[0.002 0.005 0.01 0.02 0.035 0.05 0.075]`.
- Lines 24: computes `k_rec` using `k_rec=3e8; tau_c=1e-9`.
- Lines 27: computes `pol_rf` using `pol_rf=zeros(size(fields))`.
- Lines 28: computes `pol_nz` using `pol_nz=zeros(size(fields))`.
- Lines 34: computes `pol_rf(n)` using `pol_rf(n)=geminate_pol(fields(n),k_rec,tau_c,'redfield')`.
- Lines 35: computes `pol_nz(n)` using `pol_nz(n)=geminate_pol(fields(n),k_rec,tau_c,'naka-zwan')`.

### Local helper functions

- Line 52: `geminate_pol()` — `function pol=geminate_pol(field,k_rec,tau_c,theory)`. Magnet and isotopes
  - Representative operation: `sys.magnet=field`.
  - Representative operation: `sys.isotopes={'E','E','1H'}`.

## Implementation structure

- Field dependence of geminate CIDNP from a radical pair in a viscous
- solvent, computed with Redfield theory and with the lifetime-shifted
- Nakajima-Zwanzig kernel. At a rotational correlation time of 1 ns and
- a singlet recombination rate of 1e8 Hz, the anisotropic hyperfine
- relaxation proceeds at a fair fraction of the recombination rate, and
- the pair drains before the bath decorrelates: the spectral densities
- seen by the surviving pair are lifetime-broadened, which changes the
- nuclear polarisation left in the diamagnetic product. The doubled-space
- bookkeeping follows the cidnp_geminate.m example; the low-field regime
- is motivated by the field-cycling CIDNP work of the Yurkovskaya and
- Ivanov school:
- Calculation time: minutes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pol_rf()`, `geminate_pol()`, `fields()`, `pol_nz()`, `field()`, `num2str()`, `kfigure()`, `klegend()`, `kxlabel()`, `kylabel()`, `mt2hz()`, `strcmp()`, `create()`, `basis()`, `hamiltonian()`, `assume()`.
