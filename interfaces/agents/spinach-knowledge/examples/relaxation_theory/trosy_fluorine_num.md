# examples/relaxation_theory/trosy_fluorine_num.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_fluorine_num.m`
- Signature: `trosy_fluorine_num()`
- Total lines: 102

## Purpose

Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluorine atom and its directly bonded carbon are included. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Read 3-fluorotyrosine DFT calculation; implemented by `[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'C','13C'},{'F','19F'}},[186.38 192.97])`.
- Lines 15-16: Extract coordinates and CSAs; implemented by `sys.isotopes={'19F','13C'}`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Disable startup checks; implemented by `sys.disable={'hygiene'}`.
- Lines 37-38: Magnetic field grid; implemented by `lin_freq=linspace(200,800,20)`.
- Lines 41-42: Loop over magnetic fields; implemented by `for n=1:numel(B0)`.
- Lines 44-45: Set the magnet field; implemented by `sys.magnet=B0(n)`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 54-55: States of interest; implemented by `LpF=state(spin_system,{'L+'},{1})`.
- Lines 71-72: Relaxation rates; implemented by `r2c(n)=-LpC'*R*LpC`.
- Lines 81-82: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 42: `for` loop over `n=1:numel(B0)`.

### Key state/data transformations

- Lines 12-13: computes `[~,inter_dft]` using `[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'C','13C'},{'F','19F'}},[186.38 192.97])`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'19F','13C'}`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 18: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=inter_dft.zeeman.matrix{8}`.
- Lines 19: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=inter_dft.zeeman.matrix{7}`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates=cell(2,1)`.
- Lines 21: computes `inter.coordinates{1}` using `inter.coordinates{1}=inter_dft.coordinates{8}`.
- Lines 22: computes `inter.coordinates{2}` using `inter.coordinates{2}=inter_dft.coordinates{7}`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={25e-9}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 38: computes `lin_freq` using `lin_freq=linspace(200,800,20)`.
- Lines 39: computes `B0` using `B0=2*pi*lin_freq*1e6/spin('1H')`.
- Lines 45: computes `sys.magnet` using `sys.magnet=B0(n)`.

## Implementation structure

- Transverse relaxation rate as a function of the applied magnetic
- field in a 3-fluorotyrosine labelled protein. The fluorine atom
- and its directly bonded carbon are included.
- Calculation time: minutes.
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Relaxation theory
- Basis set
- Disable startup checks
- Magnetic field grid
- Loop over magnetic fields
- Set the magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `spin()`, `create()`, `basis()`, `relaxation()`, `state()`, `r2c()`, `r2f()`, `fleft()`, `fright()`, `cleft()`, `cright()`, `kfigure()`, `kxlabel()`, `kylabel()`.
