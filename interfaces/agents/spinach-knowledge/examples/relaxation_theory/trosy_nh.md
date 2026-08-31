# examples/relaxation_theory/trosy_nh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_nh.m`
- Signature: `trosy_nh()`
- Total lines: 100

## Purpose

Transverse relaxation rate as a function of the applied magnetic field at a typical amide N-H group in a protein. Rotational cor- relation time set to 25 ns. Nitrogen CSA parameters from Nitrogen-proton bond length from DFT. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Specify coordinates and CSAs; implemented by `sys.isotopes={'1H','15N'}`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Formalism and approximation; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Disable startup checks; implemented by `sys.disable={'hygiene'}`.
- Lines 37-38: Magnetic field grid; implemented by `lin_freq=linspace(200,1500,30)`.
- Lines 41-42: Loop over magnetic fields; implemented by `for n=1:numel(B0)`.
- Lines 44-45: Set the magnet field; implemented by `sys.magnet=B0(n)`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 54-55: States of interest; implemented by `LpH=state(spin_system,{'L+'},{1})`.
- Lines 71-72: Relaxation rates; implemented by `r2c(n)=-LpN'*R*LpN`.
- Lines 81-82: Plotting; implemented by `kfigure(); plot(lin_freq',[hleft' r2h' hright'])`.

### Control flow inferred from the code

- Line 42: `for` loop over `n=1:numel(B0)`.

### Key state/data transformations

- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','15N'}`.
- Lines 17: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[ 6 0 -6]`.
- Lines 18: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[-108 62 46]`.
- Lines 19: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=pi*[0 0 0]/180`.
- Lines 20: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=pi*[0 0 -19]/180`.
- Lines 21: computes `inter.coordinates{1}` using `inter.coordinates{1}=[1.04 0.00 0.00]`.
- Lines 22: computes `inter.coordinates{2}` using `inter.coordinates{2}=[0.00 0.00 0.00]`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={25e-9}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 38: computes `lin_freq` using `lin_freq=linspace(200,1500,30)`.
- Lines 39: computes `B0` using `B0=2*pi*lin_freq*1e6/spin('1H')`.
- Lines 45: computes `sys.magnet` using `sys.magnet=B0(n)`.
- Lines 48: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Transverse relaxation rate as a function of the applied magnetic
- field at a typical amide N-H group in a protein. Rotational cor-
- relation time set to 25 ns. Nitrogen CSA parameters from
- Nitrogen-proton bond length from DFT.
- Calculation time: minutes.
- Specify coordinates and CSAs
- Relaxation theory
- Formalism and approximation
- Disable startup checks
- Magnetic field grid
- Loop over magnetic fields
- Set the magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `relaxation()`, `state()`, `r2c()`, `r2h()`, `hleft()`, `hright()`, `nleft()`, `nright()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
