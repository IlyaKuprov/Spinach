# examples/nmr_liquids/noe_four_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noe_four_spin.m`
- Signature: `noe_four_spin()`
- Total lines: 84

## Purpose

Inversion-recovery NOE effect spectrum on a simple four-spin system, with the rightmost proton signal inverted and a pulse-acquire experiment per- formed after a very long (five seconds) mixing time. Sequential NOE hops with alternating signs are clearly visible in the result. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Set the spin system; implemented by `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 21-22: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 42-43: Get thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 45-46: Start in a state with one spin inverted; implemented by `Lz1=state(spin_system,{'Lz'},{1})`.
- Lines 49-50: Evolve for five seconds under the relaxation superoperator; implemented by `rho=evolution(spin_system,1i*R,[],rho,5.0,1,'final')`.
- Lines 52-53: Subtract the non-perturbed state; implemented by `rho=rho-rho_eq`.
- Lines 55-56: Set up a pulse-acquire sequence with the resulting state set as initial; implemented by `parameters.spins={'1H'}`.
- Lines 69-70: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 72-73: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 75-76: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 78-79: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0 2.0 3.0 4.0}`.
- Lines 15: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 16: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 17: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=10`.
- Lines 18: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=10`.
- Lines 19: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]; [0.0 0.0 2.0]`.
- Lines 22: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 30: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 32: computes `inter.temperature` using `inter.temperature=298`.
- Lines 33: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `R` using `R=relaxation(spin_system)`.
- Lines 43: computes `rho_eq` using `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.

## Implementation structure

- Inversion-recovery NOE effect spectrum on a simple four-spin system, with
- the rightmost proton signal inverted and a pulse-acquire experiment per-
- formed after a very long (five seconds) mixing time. Sequential NOE hops
- with alternating signs are clearly visible in the result.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Spinach housekeeping
- Build the relaxation superoperator
- Get thermal equilibrium state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`, `evolution()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`, `scale_figure()`, `xlim()`.
