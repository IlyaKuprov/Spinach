# examples/nmr_liquids/noe_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noe_strychnine.m`
- Signature: `noe_strychnine()`
- Total lines: 84

## Purpose

Inversion-recovery NOE effect spectrum on strychnine, with the rightmost proton signal inverted and a pulse-acquire experiment performed after a 500 ms mixing time. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 15-16: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Disable Krylov propagation; implemented by `sys.disable={'krylov'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 34-35: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 44-45: Get thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 47-48: Start in a state with one spin inverted; implemented by `Lz9=state(spin_system,{'Lz'},{9})`.
- Lines 51-52: Evolve for 500 ms under the relaxation superoperator; implemented by `rho=evolution(spin_system,1i*R,[],rho,0.5,1,'final')`.
- Lines 54-55: Subtract the non-perturbed state; implemented by `rho=rho-rho_eq`.
- Lines 57-58: Set up a pulse-acquire sequence with the resulting state set as initial; implemented by `parameters.spins={'1H'}`.
- Lines 71-72: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 74-75: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 77-78: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 80-81: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 24: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 25: computes `bas.space_level` using `bas.space_level=3`.
- Lines 28: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 30: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 31: computes `inter.temperature` using `inter.temperature=298`.
- Lines 32: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 35: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `R` using `R=relaxation(spin_system)`.
- Lines 45: computes `rho_eq` using `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 48: computes `Lz9` using `Lz9=state(spin_system,{'Lz'},{9})`.
- Lines 49: computes `rho` using `rho=rho_eq-2*Lz9*(Lz9'*rho_eq)/norm(Lz9)^2`.

## Implementation structure

- Inversion-recovery NOE effect spectrum on strychnine, with the rightmost
- proton signal inverted and a pulse-acquire experiment performed after a
- 500 ms mixing time.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Disable Krylov propagation
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Build the relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`, `evolution()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
