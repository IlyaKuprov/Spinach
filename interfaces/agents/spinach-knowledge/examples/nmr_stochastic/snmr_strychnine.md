# examples/nmr_stochastic/snmr_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_stochastic/snmr_strychnine.m`
- Signature: `snmr_strychnine()`
- Total lines: 107

## Purpose

A Primas-style stochastic NMR experiment on strychnine. The calculation requires an NVidia Titan V GPU at a minimum. Calculation time: minutes

## Physical / mathematical content

- Stochastic NMR examples. These scripts model random processes, trajectories, or stochastic Liouville dynamics and connect fluctuating Hamiltonians or transport processes to ensemble-averaged observables.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Algorithmic options; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 32-36: Use GPU arithmetic sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Get the Hamiltonian; implemented by `H0=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 42-43: Get the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 45-46: Stochastic process parameters; implemented by `dt=1e-5`.
- Lines 50-51: Control operators; implemented by `Hx=operator(spin_system,'Lx','1H')`.
- Lines 54-55: Uniformly distributed noise; implemented by `Cx=omega_max*(2*rand(nsteps,1)-1)`.
- Lines 58-59: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 61-62: GPU uploads; implemented by `H0=gpuArray(H0); R=gpuArray(R)`.
- Lines 65-66: Trajectory calculation; implemented by `report(spin_system,'Computing trajectory ')`.
- Lines 71-72: Evolution generator; implemented by `G=H0+1i*R+Cx(n)*Hx+Cy(n)*Hy`.
- Lines 74-75: Time step; implemented by `traj(:,n+1)=step(spin_system,G,traj(:,n),dt)`.
- Lines 79-80: Performance report; implemented by `disp(['Steps per second: ' num2str(nsteps/toc)])`.

### Control flow inferred from the code

- Line 69: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 28: computes `inter.temperature` using `inter.temperature=298`.
- Lines 29: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 30: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `H0` using `H0=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 43: computes `R` using `R=relaxation(spin_system)`.
- Lines 46: computes `dt` using `dt=1e-5`.
- Lines 47: computes `omega_max` using `omega_max=1e2`.
- Lines 48: computes `nsteps` using `nsteps=10000`.

## Implementation structure

- A Primas-style stochastic NMR experiment on strychnine. The calculation
- requires an NVidia Titan V GPU at a minimum.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Relaxation theory parameters
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Get the Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `operator()`, `equilibrium()`, `gpuArray()`, `report()`, `traj()`, `step()`, `num2str()`, `gather()`, `state()`, `kfigure()`.
