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
