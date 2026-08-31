# examples/relaxation_theory/from_md/ngce_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/from_md/ngce_test.m`
- Signature: `ngce_test()`
- Total lines: 77

## Purpose

Test of the numerical integral route to the Redfield relaxation superoperator against the analytical results for isotropic rota- tional diffusion. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=9.4`.
- Lines 18-19: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Get Redfield relaxation matrix; implemented by `R_red=real(relaxation(spin_system))`.
- Lines 35-36: Get lab frame Hamiltonian components; implemented by `[H0,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 38-39: Get a random walk on a sphere; implemented by `tau_c=inter.tau_c{1}; dt=tau_c/25`.
- Lines 42-43: Get Hamiltonian trajectory; implemented by `H1=cell(size(eulers,1),1)`.
- Lines 48-49: Compute the relaxation matrix using GCE; implemented by `[R_gce,dR_gce]=ngce(spin_system,H0,H1,dt,tau_c,1e-3)`.
- Lines 51-52: Get the answers to the user; implemented by `disp('Relaxation superoperator, analytical'); disp(full(R_red))`.
- Lines 56-57: Get relaxation rates; implemented by `gce_rates_stdm=diag(dR_gce)`.
- Lines 61-62: Draw a diagonal straight line; implemented by `min_rate=min([gce_rates; red_rates])`.
- Lines 68-71: Plot 95% confidence bands for the rates; implemented by `errorbar(gce_rates,red_rates,2*gce_rates_stdm, 'horizontal','LineStyle','none', 'Marker','o','Color','r'); hold on`.

### Control flow inferred from the code

- Line 44: `parfor` loop over `n=1:size(eulers,1)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 15: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 19: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 20: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 22: computes `inter.tau_c` using `inter.tau_c={1.0e-10}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `R_red` using `R_red=real(relaxation(spin_system))`.
- Lines 36: computes `[H0,Q]` using `[H0,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 39: computes `tau_c` using `tau_c=inter.tau_c{1}; dt=tau_c/25`.
- Lines 40: computes `eulers` using `eulers=rwalk(100000,tau_c,dt)`.
- Lines 43: computes `H1` using `H1=cell(size(eulers,1),1)`.
- Lines 45: computes `H1{n}` using `H1{n}=orientation(Q,eulers(n,:))`.
- Lines 49: computes `[R_gce,dR_gce]` using `[R_gce,dR_gce]=ngce(spin_system,H0,H1,dt,tau_c,1e-3)`.
- Lines 57: computes `gce_rates_stdm` using `gce_rates_stdm=diag(dR_gce)`.

## Implementation structure

- Test of the numerical integral route to the Redfield relaxation
- superoperator against the analytical results for isotropic rota-
- tional diffusion.
- Calculation time: minutes.
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Get Redfield relaxation matrix
- Get lab frame Hamiltonian components
- Get a random walk on a sphere
- Get Hamiltonian trajectory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `hamiltonian()`, `assume()`, `rwalk()`, `orientation()`, `eulers()`, `ngce()`, `STD()`, `kfigure()`, `errorbar()`, `kxlabel()`, `kylabel()`.
