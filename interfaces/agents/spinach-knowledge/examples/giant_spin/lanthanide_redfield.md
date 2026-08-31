# examples/giant_spin/lanthanide_redfield.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/lanthanide_redfield.m`
- Signature: `lanthanide_redfield()`
- Total lines: 63

## Purpose

Relaxation rate of Gd(III) as a function of zero-field splitting, computed using Redfield theory. The correla- tion time (1 fs) refers to vibrational dynamics of the ligand cage. Calculation time: minutes

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Spin system properties; implemented by `sys.isotopes={'E8'}`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=9.40`.
- Lines 20-21: Relaxation parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Linearly spaced B20 in cm^-1; implemented by `b20=linspace(0.1,10,10)'`.
- Lines 33-34: Loop over B20 values; implemented by `for n=1:numel(b20)`.
- Lines 36-37: Giant spin Hamiltonian parameters; implemented by `inter.giant.coeff={{[0 0 0],[0 0 icm2hz(b20(n)) 0 0]}}`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 47-48: Longitudinal and tranverse relaxation rates; implemented by `Lp=state(spin_system,'L+','E8')`.
- Lines 55-56: Do the plotting; implemented by `kfigure(); plot(b20,[T1; T2],'o')`.

### Control flow inferred from the code

- Line 34: `for` loop over `n=1:numel(b20)`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 18: computes `sys.magnet` using `sys.magnet=9.40`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.tau_c` using `inter.tau_c={1e-15}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `b20` using `b20=linspace(0.1,10,10)'`.
- Lines 37: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],[0 0 icm2hz(b20(n)) 0 0]}}`.
- Lines 38: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0]}}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `R` using `R=relaxation(spin_system)`.
- Lines 48: computes `Lp` using `Lp=state(spin_system,'L+','E8')`.
- Lines 49: computes `Lz` using `Lz=state(spin_system,'Lz','E8')`.
- Lines 50: computes `T1(n)` using `T1(n)=-1/((Lz'*R*Lz)/(Lz'*Lz))`.
- Lines 51: computes `T2(n)` using `T2(n)=-1/((Lp'*R*Lp)/(Lp'*Lp))`.

## Implementation structure

- Relaxation rate of Gd(III) as a function of zero-field
- splitting, computed using Redfield theory. The correla-
- tion time (1 fs) refers to vibrational dynamics of the
- ligand cage.
- Calculation time: minutes
- Spin system properties
- Magnet field
- Relaxation parameters
- Basis set
- Linearly spaced B20 in cm^-1
- Loop over B20 values
- Giant spin Hamiltonian parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `b20()`, `create()`, `basis()`, `relaxation()`, `state()`, `kfigure()`, `set()`, `kxlabel()`, `kylabel()`.
