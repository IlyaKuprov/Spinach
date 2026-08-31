# examples/optimal_control/magic_pulse_phase.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/magic_pulse_phase.m`
- Signature: `magic_pulse_phase()`
- Total lines: 135

## Purpose

A template file for the "magic pulse" optimisations. The term refers to a family of broadband NMR pulses that are tolerant to resonance offsets and power calibration errors: Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz) to be negligible. The latter requ

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Set the magnetic field; implemented by `sys.magnet=28.18`.
- Lines 25-27: Put 100 non-interacting spins at equal intervals within the [-100,+100] ppm chemical shift range; implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 33-35: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Set up spin states; implemented by `Sx=state(spin_system,'Lx','13C')`.
- Lines 52-53: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 56-57: Get the drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 59-60: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 72-74: Plotting options; implemented by `control.plotting={'phi_controls','xy_controls', 'robustness','spectrogram'}`.
- Lines 76-77: Initial guess; implemented by `guess=(pi/5)*randn(1,60)`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 82-83: Run LBFGS GRAPE pulse optimisation; implemented by `phi_profile=fmaxnewton(spin_system,@grape_phase,guess)`.
- Lines 85-86: Get Cartesian components of the pulse; implemented by `amp_profile=mean(control.pwr_levels)*control.amplitudes`.
- Lines 89-90: Simulate the optimised pulse; implemented by `rho_init=state(spin_system,'Lz','13C')`.
- Lines 94-95: Set acquisition parameters; implemented by `parameters.spins={'13C'}`.
- Lines 106-107: Simulate the free induction decay; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 109-110: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.
- Lines 112-113: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.

### Control flow inferred from the code

- Line 28: `for` loop over `n=1:n_spins`.

### Key state/data transformations

- Lines 23: computes `sys.magnet` using `sys.magnet=28.18`.
- Lines 27: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 29: computes `sys.isotopes{n}` using `sys.isotopes{n}='13C'`.
- Lines 31: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins))`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 37: computes `bas.space_level` using `bas.space_level=1`.
- Lines 38: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `Sx` using `Sx=state(spin_system,'Lx','13C')`.
- Lines 46: computes `Sy` using `Sy=state(spin_system,'Ly','13C')`.
- Lines 47: computes `Sz` using `Sz=state(spin_system,'Lz','13C')`.
- Lines 53: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 54: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 57: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 60: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 61: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 62: computes `control.drifts` using `control.drifts={{H}}`.

## Implementation structure

- A template file for the "magic pulse" optimisations. The term refers to
- a family of broadband NMR pulses that are tolerant to resonance offsets
- and power calibration errors:
- Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The
- pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and
- must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz)
- to be negligible. The latter requirement caps the duration at 1/100*J
- = 50 us. The pulse must accomplish the following transfers: {Lz -> Lx,
- Ly -> Ly, Lx -> -Lz}. A realistically achievable nutation frequency is
- between 50 kHz and 70 kHz across the RF coil.
- Calculation time: minutes.
- Set the magnetic field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `polar2cartesian()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `figure()`, `subplot()`.
