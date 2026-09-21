# examples/giant_spin/case_studies/mn3_trimer_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/mn3_trimer_magn.m`
- Signature: `mn3_trimer_magn()`
- Total lines: 88

## Purpose

Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a linear trimer of three S=5/2 manganese ions with isotropic exchange between neighbours and an axial plus rhombic zero-field splitting on every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The full 216-state Hilbert space is used; the paper solves the same problem in 

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Magnet must be 1 Tesla, the field is set by the sweep; implemented by `sys.magnet=1.0`.
- Lines 22-23: Parallel pool size; implemented by `sys.parallel={'processes',4}`.
- Lines 25-26: Three S=5/2 spins with g=2, phonon bath at 0.6 K; implemented by `sys.isotopes={'E6','E6','E6'}`.
- Lines 30-31: Isotropic exchange, J=-2.42 cm^-1 in the H=-2*J*S1*S2 convention of the paper; implemented by `inter.coupling.matrix=cell(3,3)`.
- Lines 35-36: Zero-field splitting, D=0.167 cm^-1 and E=0.040 cm^-1 on every ion; implemented by `for n=1:3`.
- Lines 40-41: Formalism and basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Total S_z and the spin-phonon coupling operator between adjacent total S_z states; implemented by `Sz=full(operator(spin_system,'Lz','E6')); msz=diag(Sz)`.
- Lines 52-53: Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units; implemented by `parameters.phonon_alpha=2`.
- Lines 56-57: Observable: total magnetic moment along Z in Bohr magnetons; implemented by `parameters.coil=-2.0*Sz`.
- Lines 59-60: Sweep: 50 T/ms to 10 T in 20 ns stairs, output every 100 stairs; implemented by `parameters.field_prof=@(t) 5e4*t`.
- Lines 63-64: Single crystal in the frame of the zero-field splitting tensors; implemented by `parameters.spins={'E6'}; parameters.orientation=[0 0 0]`.
- Lines 67-68: Run the simulation; implemented by `answer=crystal(spin_system,@pulsed_field,parameters,'labframe')`.
- Lines 70-71: Thermal equilibrium magnetisation at the same fields; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0])`.
- Lines 79-80: Plot the sweep and the equilibrium curves; implemented by `kfigure(); plot(answer.field,answer.obs); hold on; plot(answer.field,m_eq,'--'); hold off`.
- Lines 84-85: Save the curves; implemented by `save('mn3_trimer_magn.mat','answer','m_eq')`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:3`.
- Line 73: `for` loop over `k=1:numel(m_eq)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 23: computes `sys.parallel` using `sys.parallel={'processes',4}`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'E6','E6','E6'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0 2.0 2.0}`.
- Lines 28: computes `inter.temperature` using `inter.temperature=0.6`.
- Lines 31: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(3,3)`.
- Lines 32: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=-2*icm2hz(-2.42)*eye(3)`.
- Lines 33: computes `inter.coupling.matrix{2,3}` using `inter.coupling.matrix{2,3}=-2*icm2hz(-2.42)*eye(3)`.
- Lines 37: computes `inter.coupling.matrix{n,n}` using `inter.coupling.matrix{n,n}=zfs2mat(icm2hz(0.167),icm2hz(0.040),0,0,0)`.
- Lines 41: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `Sz` using `Sz=full(operator(spin_system,'Lz','E6')); msz=diag(Sz)`.
- Lines 53: computes `parameters.phonon_alpha` using `parameters.phonon_alpha=2`.
- Lines 54: computes `parameters.phonon_i0` using `parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2`.
- Lines 57: computes `parameters.coil` using `parameters.coil=-2.0*Sz`.
- Lines 60: computes `parameters.field_prof` using `parameters.field_prof=@(t) 5e4*t`.
- Lines 61: computes `parameters.timestep` using `parameters.timestep=2e-8; parameters.nsteps=1e4; parameters.nout=100`.

## Implementation structure

- Pulsed-field magnetisation of the (CH6N3)2MnCl4 molecular crystal, a
- linear trimer of three S=5/2 manganese ions with isotropic exchange
- between neighbours and an axial plus rhombic zero-field splitting on
- every ion, at 0.6 K under a 50 T/ms sweep to 10 T with spin-phonon
- relaxation in the generalised Lindblad form of Saito and Miyashita.
- The full 216-state Hilbert space is used; the paper solves the same
- problem in 16-state and 26-state effective bases. The thermal equi-
- librium magnetisation is plotted for comparison. Reproduces Figure 6
- Calculation time: minutes
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size
- Three S=5/2 spins with g=2, phonon bath at 0.6 K

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `operator()`, `double()`, `crystal()`, `hamiltonian()`, `assume()`, `orientation()`, `m_eq()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `save()`.
