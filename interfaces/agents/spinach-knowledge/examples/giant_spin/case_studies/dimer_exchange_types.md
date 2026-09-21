# examples/giant_spin/case_studies/dimer_exchange_types.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/dimer_exchange_types.m`
- Signature: `dimer_exchange_types()`
- Total lines: 96

## Purpose

Pulsed-field magnetisation of a dimer of two S=1/2 spins with four types of exchange coupling tensor: isotropic, two anisotropic, and antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin- phonon relaxation in the generalised Lindblad form of Saito and Miyashita. The out-of-equilibrium curves are compared with the thermal equilibrium magnetisation. Reproduces Figure 4 of Calculation time: minutes

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet must be 1 Tesla, the field is set by the sweep; implemented by `sys.magnet=1.0`.
- Lines 19-20: Parallel pool size; implemented by `sys.parallel={'processes',4}`.
- Lines 22-23: Two electron spins with g=2; implemented by `sys.isotopes={'E','E'}`.
- Lines 26-27: Temperature of the phonon bath; implemented by `inter.temperature=0.2`.
- Lines 29-30: Formalism and basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 33-37: Exchange coupling tensors of the paper, cm^-1, in its H=-2*S1*J*S2 convention; implemented by `tensors={[0.2 0 0; 0 0.2 0; 0 0 0.2], [0.2 0 0; 0 0 0; 0 0 0 ], [0 0 0; 0 0 0; 0 0 0.2], [0 0.2 0.2; -0.2 0 0.2; -0.2 -0.2 0]}`.
- Lines 40-41: Spin-phonon coupling operator: unit elements between adjacent total S_z states; implemented by `Sz=kron(full(stevens(2,1,0)),eye(2))+kron(eye(2),full(stevens(2,1,0))); msz=diag(Sz)`.
- Lines 44-45: Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-10 ps/rad) in rad/s units; implemented by `parameters.phonon_alpha=2`.
- Lines 48-49: Observable: total magnetic moment along Z in Bohr magnetons; implemented by `parameters.coil=-2.0*Sz`.
- Lines 51-52: Sweep: 10 T/ms to 1 T in 10 ns stairs, output every 10 stairs; implemented by `parameters.field_prof=@(t) 1e4*t`.
- Lines 55-56: Single crystal in the frame of the exchange tensor; implemented by `parameters.spins={'E'}; parameters.orientation=[0 0 0]`.
- Lines 59-60: Loop over the exchange tensors; implemented by `kfigure(); scale_figure([2.0 1.6]); answers=cell(1,4)`.
- Lines 63-64: Spinach coupling convention is S1*A*S2 with A in Hz; implemented by `inter.coupling.matrix=cell(2,2)`.
- Lines 67-68: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 71-72: Run the simulation; implemented by `answers{n}=crystal(spin_system,@pulsed_field,parameters,'labframe')`.
- Lines 74-75: Thermal equilibrium magnetisation at the same fields; implemented by `[I,Q]=hamiltonian(assume(spin_system,'labframe')); H0=I+orientation(Q,[0 0 0])`.
- Lines 84-85: Plot the sweep and the equilibrium curves; implemented by `subplot(2,2,n); plot(answers{n}.field,answers{n}.obs); hold on`.
- Lines 92-93: Save the curves; implemented by `save('dimer_exchange_types.mat','answers','labels')`.

### Control flow inferred from the code

- Line 61: `for` loop over `n=1:4`.
- Line 77: `for` loop over `k=1:numel(m_eq)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 20: computes `sys.parallel` using `sys.parallel={'processes',4}`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0 2.0}`.
- Lines 27: computes `inter.temperature` using `inter.temperature=0.2`.
- Lines 30: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34-37: computes `tensors` using `tensors={[0.2 0 0; 0 0.2 0; 0 0 0.2], [0.2 0 0; 0 0 0; 0 0 0 ], [0 0 0; 0 0 0; 0 0 0.2], [0 0.2 0.2; -0.2 0 0.2; -0.2 -0.2 0]}`.
- Lines 38: computes `labels` using `labels={'isotropic','anisotropic, $J_{xx}$','anisotropic, $J_{zz}$','antisymmetric'}`.
- Lines 41: computes `Sz` using `Sz=kron(full(stevens(2,1,0)),eye(2))+kron(eye(2),full(stevens(2,1,0))); msz=diag(Sz)`.
- Lines 45: computes `parameters.phonon_alpha` using `parameters.phonon_alpha=2`.
- Lines 46: computes `parameters.phonon_i0` using `parameters.phonon_i0=1e2*1e-10*1e12*(1e-12)^2*0.1883651568463003^2`.
- Lines 49: computes `parameters.coil` using `parameters.coil=-2.0*Sz`.
- Lines 52: computes `parameters.field_prof` using `parameters.field_prof=@(t) 1e4*t`.
- Lines 53: computes `parameters.timestep` using `parameters.timestep=1e-8; parameters.nsteps=1e4; parameters.nout=10`.
- Lines 56: computes `parameters.spins` using `parameters.spins={'E'}; parameters.orientation=[0 0 0]`.
- Lines 57: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 60: computes `kfigure(); scale_figure([2.0 1.6]); answers` using `kfigure(); scale_figure([2.0 1.6]); answers=cell(1,4)`.

## Implementation structure

- Pulsed-field magnetisation of a dimer of two S=1/2 spins with four
- types of exchange coupling tensor: isotropic, two anisotropic, and
- antisymmetric, at 0.2 K under a 10 T/ms sweep to 1 T, with spin-
- phonon relaxation in the generalised Lindblad form of Saito and
- Miyashita. The out-of-equilibrium curves are compared with the
- thermal equilibrium magnetisation. Reproduces Figure 4 of
- Calculation time: minutes
- Magnet must be 1 Tesla, the field is set by the sweep
- Parallel pool size
- Two electron spins with g=2
- Temperature of the phonon bath
- Formalism and basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `stevens()`, `double()`, `kfigure()`, `scale_figure()`, `icm2hz()`, `create()`, `basis()`, `crystal()`, `hamiltonian()`, `assume()`, `orientation()`, `m_eq()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
