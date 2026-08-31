# examples/relaxation_theory/inv_rec_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/inv_rec_2.m`
- Signature: `inv_rec_2()`
- Total lines: 97

## Purpose

An example of inversion recovery experiment simulation for a strychnine spin system. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Disable Krylov propagation; implemented by `sys.disable={'krylov'}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 32-33: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Aquisition parameters; implemented by `parameters.sweep=6500`.
- Lines 47-48: Set up different recovery delays; implemented by `mixing_time=[0.01; 0.1; 0.5; 1; 5; 10]`.
- Lines 50-51: Isotropic thermal equilibrium; implemented by `parameters.rho_eq=equilibrium(spin_system)`.
- Lines 53-54: Detection state; implemented by `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 56-58: Rotating frame Liouvillian; implemented by `L=hamiltonian(assume(spin_system,'nmr'))+ 1i*relaxation(spin_system)`.
- Lines 60-61: Pulse operator; implemented by `Ly=operator(spin_system,'Ly','1H')`.
- Lines 63-64: Get figure going; implemented by `kfigure(); scale_figure([2.5 1.8])`.
- Lines 66-67: Set up loop; implemented by `for n=1:numel(mixing_time)`.
- Lines 69-70: Run different recovery delays; implemented by `parameters.tau=mixing_time(n)`.
- Lines 72-73: Apply the 180 degree inversion pulse; implemented by `rho=step(spin_system,Ly,parameters.rho_eq,pi)`.
- Lines 75-76: Evolution; implemented by `rho=evolution(spin_system,L,[],rho,parameters.tau,1,'final')`.

### Control flow inferred from the code

- Line 67: `for` loop over `n=1:numel(mixing_time)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=3`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 29: computes `inter.temperature` using `inter.temperature=298`.
- Lines 30: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 33: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=6500`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=8192`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=65536`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- An example of inversion recovery experiment simulation
- for a strychnine spin system.
- Calculation time: minutes.
- Read the spin system properties
- Magnet field
- Disable Krylov propagation
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Aquisition parameters
- Set up different recovery delays

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `equilibrium()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `operator()`, `kfigure()`, `scale_figure()`, `mixing_time()`, `step()`, `evolution()`, `apodisation()`, `fftshift()`.
