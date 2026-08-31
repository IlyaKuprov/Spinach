# examples/optimal_control/static_powder_control.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/static_powder_control.m`
- Signature: `static_powder_control()`
- Total lines: 138

## Purpose

Optimal control optimisation for a pulse that is designed to set deuterium magnetisation in a -CD3 group of alanine up for perfect rephasing 100 microseconds after the pulse is finished. The system is a powder (100 orientations) with a B1 dist- ribution (from 46 to 54 kHz per channel) and transmitter offset error within 1 kHz of the chemical shift. Goodwin's very efficient version of the GRAPE Hessian al- gorithm is 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: 600 MHz magnet; implemented by `sys.magnet=14.1`.
- Lines 25-26: Isotopes; implemented by `sys.isotopes={'2H'}`.
- Lines 28-29: Alanine CD3 NQI parameters; implemented by `inter.coupling.matrix{1,1}=anas2mat(0,40e3,0,0,0,0)`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Powder context parameters; implemented by `parameters.spins={'2H'}`.
- Lines 47-48: Drift Liouvillians for the ensemble; implemented by `control.drifts=drifts(spin_system,@powder,parameters,'labframe')`.
- Lines 50-51: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,'Lz','2H')`.
- Lines 54-55: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,'Lx','2H')`.
- Lines 58-59: Get deuterium control operators; implemented by `Dx=operator(spin_system,'Lx','2H')`.
- Lines 62-63: Get deuterium offset operator; implemented by `Dz=operator(spin_system,'Lz','2H')`.
- Lines 65-66: Define control parameters; implemented by `control.isotopes={'2H'}`.
- Lines 81-82: Plotting options; implemented by `control.plotting={'xy_controls','robustness','spectrogram'}`.
- Lines 84-85: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 87-88: Waveform guess; implemented by `guess=randn(2,100)/10`.
- Lines 90-91: Run the optimisation; implemented by `pulse_profile=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 93-94: Get Cartesian components of the pulse; implemented by `CLx=mean(control.pwr_levels)*pulse_profile(1,:)`.
- Lines 97-98: Run a test calculation; implemented by `fid_optim=zeros([500 1],'like',1i)`.

### Control flow inferred from the code

- Line 100: `parfor` loop over `n=1:numel(control.drifts)`.

### Key state/data transformations

- Lines 23: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'2H'}`.
- Lines 29: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=anas2mat(0,40e3,0,0,0,0)`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 41: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 42: computes `parameters.offset` using `parameters.offset=0`.
- Lines 43: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 44: computes `parameters.rframes` using `parameters.rframes={{'2H',2}}`.
- Lines 45: computes `parameters.verbose` using `parameters.verbose=0`.
- Lines 48: computes `control.drifts` using `control.drifts=drifts(spin_system,@powder,parameters,'labframe')`.
- Lines 51: computes `rho_init` using `rho_init=state(spin_system,'Lz','2H')`.
- Lines 55: computes `rho_targ` using `rho_targ=state(spin_system,'Lx','2H')`.
- Lines 59: computes `Dx` using `Dx=operator(spin_system,'Lx','2H')`.
- Lines 60: computes `Dy` using `Dy=operator(spin_system,'Ly','2H')`.
- Lines 63: computes `Dz` using `Dz=operator(spin_system,'Lz','2H')`.

## Implementation structure

- Optimal control optimisation for a pulse that is designed
- to set deuterium magnetisation in a -CD3 group of alanine
- up for perfect rephasing 100 microseconds after the pulse
- is finished.
- The system is a powder (100 orientations) with a B1 dist-
- ribution (from 46 to 54 kHz per channel) and transmitter
- offset error within 1 kHz of the chemical shift.
- Goodwin's very efficient version of the GRAPE Hessian al-
- gorithm is used because propagator dimensions are small;
- it yields a sophisticated kind of spin echo.
- Calculation time: minutes
- 600 MHz magnet

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `anas2mat()`, `create()`, `basis()`, `drifts()`, `state()`, `operator()`, `optimcon()`, `fmaxnewton()`, `pulse_profile()`, `shaped_pulse_xy()`, `evolution()`, `kfigure()`, `subplot()`, `kxlabel()`, `xlim()`, `fid_optim()`.
