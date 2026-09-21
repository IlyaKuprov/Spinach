# examples/giant_spin/case_studies/ho_pzdo4_profiles.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_profiles.m`
- Signature: `ho_pzdo4_profiles()`
- Total lines: 117

## Purpose

Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under four magnetic field profiles: linear sweep, piecewise linear sweep, monotone cubic spline through a measured 65 T short pulse, and a sinu- soidal field at the clock transition frequency. Spin-phonon relaxation is the generalised Lindblad dissipator of Saito and Miyashita with a s

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Crystal field parameters, cm^-1, ranks 2 to 12 in Stevens operator convention; implemented by `[ks,qs,bkq]=ho_pzdo4_params()`.
- Lines 25-26: Convert Stevens coefficients into spherical tensor coefficients, Hz, rank by rank; implemented by `coeff=cell(1,12); euler=cell(1,12)`.
- Lines 32-33: Magnet must be 1 Tesla, the field is set by the sweep; implemented by `sys.magnet=1.0`.
- Lines 35-36: Parallel pool size; implemented by `sys.parallel={'processes',4}`.
- Lines 38-39: J=8 giant spin, effective g-factor 1.24; implemented by `sys.isotopes={'E17'}`.
- Lines 44-45: Formalism and basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 48-49: Spin-phonon coupling operator: unit elements between adjacent m_J states; implemented by `Jz=full(stevens(17,1,0)); mj=diag(Jz)`.
- Lines 52-53: Super-Ohmic bath, lambda^2*I0 of the paper (lambda=10 cm^-1, I0=1e-14 ps/rad) in rad/s units; implemented by `parameters.phonon_alpha=2`.
- Lines 56-57: Observable: magnetic moment along Z in Bohr magnetons; implemented by `parameters.coil=-1.24*Jz`.
- Lines 59-60: Single crystal, crystal field frame aligned with the laboratory frame; implemented by `parameters.spins={'E17'}; parameters.orientation=[0 0 0]`.
- Lines 63-67: Measured 65 T pulse of the paper, ms and Tesla, 76 points, for the monotone spline profile; implemented by `pulse_t=[0 0.2 0.3 0.4 0.5 0.7 0.8 1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.2 3.3 3.4 3.5 3.7 3.8 4 4.1 4.2 4.3 4.4 4.7 4.8 5 5.1 5.2…`.
- Lines 79-83: Four field profiles, Tesla as a function of time in seconds; implemented by `profiles={@(t) 1e4*t, @(t) interp1([0 1e-6 1e-5 1e-4 1e-3],[0 0.1 1 5 10],t,'linear'), @(t) pchip(pulse_t,pulse_b,t), @(t) 0.1*sin(0.134124264765e12*t)}`.
- Lines 88-89: Loop over the profiles; implemented by `kfigure(); scale_figure([2.0 1.6]); answers=cell(1,4)`.
- Lines 92-93: Spinach housekeeping at the temperature of the panel; implemented by `inter.temperature=temps(n)`.
- Lines 97-98: Sweep parameters of the panel; implemented by `parameters.field_prof=profiles{n}`.
- Lines 103-104: Run the simulation; implemented by `answers{n}=crystal(spin_system,@pulsed_field,parameters,'labframe')`.
- Lines 106-107: Plot the field and the magnetisation against time; implemented by `subplot(2,2,n); yyaxis left; plot(answers{n}.t*tscale(n),answers{n}.field); kylabel('Field, Tesla')`.
- Lines 113-114: Save the curves; implemented by `save('ho_pzdo4_profiles.mat','answers','labels')`.

### Control flow inferred from the code

- Line 27: `for` loop over `k=1:12`.
- Line 90: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 23: computes `[ks,qs,bkq]` using `[ks,qs,bkq]=ho_pzdo4_params()`.
- Lines 26: computes `coeff` using `coeff=cell(1,12); euler=cell(1,12)`.
- Lines 29: computes `coeff{k}` using `coeff{k}=stev2sph(k,icm2hz(stev)); euler{k}=[0 0 0]`.
- Lines 33: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 36: computes `sys.parallel` using `sys.parallel={'processes',4}`.
- Lines 39: computes `sys.isotopes` using `sys.isotopes={'E17'}`.
- Lines 40: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.24}`.
- Lines 41: computes `inter.giant.coeff` using `inter.giant.coeff={coeff}`.
- Lines 42: computes `inter.giant.euler` using `inter.giant.euler={euler}`.
- Lines 45: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 46: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 49: computes `Jz` using `Jz=full(stevens(17,1,0)); mj=diag(Jz)`.
- Lines 53: computes `parameters.phonon_alpha` using `parameters.phonon_alpha=2`.
- Lines 54: computes `parameters.phonon_i0` using `parameters.phonon_i0=1e2*1e-14*1e12*(1e-12)^2*0.1883651568463003^2`.
- Lines 57: computes `parameters.coil` using `parameters.coil=-1.24*Jz`.
- Lines 60: computes `parameters.spins` using `parameters.spins={'E17'}; parameters.orientation=[0 0 0]`.
- Lines 61: computes `parameters.needs` using `parameters.needs={'zeeman_op'}`.
- Lines 64-67: computes `pulse_t` using `pulse_t=[0 0.2 0.3 0.4 0.5 0.7 0.8 1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.2 3.3 3.4 3.5 3.7 3.8 4 4.1 4.2 4.3 4.4 4.7 4.8 5 5.1 5.2…`.

## Implementation structure

- Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework,
- a J=8 giant spin with a crystal field to twelfth spherical rank, under
- four magnetic field profiles: linear sweep, piecewise linear sweep,
- monotone cubic spline through a measured 65 T short pulse, and a sinu-
- soidal field at the clock transition frequency. Spin-phonon relaxation
- is the generalised Lindblad dissipator of Saito and Miyashita with a
- super-Ohmic phonon bath. Reproduces Figure 2 of
- with the crystal field parameters, g-factor, temperatures, spectral
- density, sweep profiles, and stair widths of that paper. As in the
- paper, the first three profiles are propagated for 1 ms (the first
- millisecond of the measured 10 ms pulse), the sinusoid for 140 ps.
- Calculation time: minutes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ho_pzdo4_params()`, `stev()`, `bkq()`, `stev2sph()`, `icm2hz()`, `stevens()`, `double()`, `pchip()`, `kfigure()`, `scale_figure()`, `temps()`, `create()`, `basis()`, `steps()`, `nsteps()`, `nout()`.
