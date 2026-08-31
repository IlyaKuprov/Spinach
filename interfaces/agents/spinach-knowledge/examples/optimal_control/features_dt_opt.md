# examples/optimal_control/features_dt_opt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_dt_opt.m`
- Signature: `features_dt_opt()`
- Total lines: 110

## Purpose

Optimisation of slice durations in a composite inversion pulse with specified amplitudes, phases, and a constrain- ed overall duration. The initial guess is 270(-x)360(x)90(y)270(-y)360(y)90(x) [Fig. 3] from https://doi.org/10.1016/0022-2364(83)90133-6 --the optimisation demonstrates that a slightly better pulse of the same power and duration exists. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Set the magnetic field; implemented by `sys.magnet=14.1`.
- Lines 20-22: Put 100 non-interacting spins at equal intervals over the area that needs to be affected by the pulse (25 kHz either side); implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 28-30: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Set up initial and target states; implemented by `rho_init=state(spin_system,'Lz','13C')`.
- Lines 44-45: Get the control operators; implemented by `controls{1}=operator(spin_system,'Lx','13C')`.
- Lines 48-49: Get the drift Hamiltonian; implemented by `drift=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 51-52: Specify the pulse sequence (rad/s); implemented by `waveform=2*pi*2.5e4*[-1 1 0 0 0 1`.
- Lines 56-58: Optimisation options; implemented by `options=optimoptions('fmincon','SpecifyObjectiveGradient',true, 'Display','iter','HessianApproximation','lbfgs')`.
- Lines 60-62: Use Matlab's optimiser with a total duration constraint; implemented by `target_fun=@(x)tgrape(spin_system,drift,controls,waveform, x,time_unit,rho_init,rho_targ)`.
- Lines 66-67: Apply the old and the new pulse; implemented by `rho_old=rho_init; rho_new=rho_init`.
- Lines 76-77: Set acquisition parameters; implemented by `parameters.spins={'13C'}`.
- Lines 87-88: Run the simulation for the initial guess; implemented by `parameters.rho0=step(spin_system,Ly,rho_old,pi/2)`.
- Lines 93-94: Run the simulation for the optimised pulse; implemented by `parameters.rho0=step(spin_system,Ly,rho_new,pi/2)`.
- Lines 99-100: Plot the spectrum before and after optimisation; implemented by `kfigure(); scale_figure([2.0 0.75])`.
- Lines 106-107: Write the old and the new pulses to the console; implemented by `disp('Old and new durations (microseconds):'); disp([dt_old dt_new])`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:n_spins`.
- Line 69: `for` loop over `n=1:numel(dt_new)`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 22: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 24: computes `sys.isotopes{n}` using `sys.isotopes{n}='13C'`.
- Lines 26: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-166,166,n_spins))`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 32: computes `bas.space_level` using `bas.space_level=1`.
- Lines 33: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `rho_init` using `rho_init=state(spin_system,'Lz','13C')`.
- Lines 42: computes `rho_targ` using `rho_targ=rho_init`.
- Lines 45: computes `controls{1}` using `controls{1}=operator(spin_system,'Lx','13C')`.
- Lines 46: computes `controls{2}` using `controls{2}=operator(spin_system,'Ly','13C')`.
- Lines 49: computes `drift` using `drift=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 52: computes `waveform` using `waveform=2*pi*2.5e4*[-1 1 0 0 0 1`.
- Lines 54: computes `dt_old` using `dt_old=10*[3; 4; 1; 3; 4; 1]; time_unit=1e-6`.
- Lines 57-58: computes `options` using `options=optimoptions('fmincon','SpecifyObjectiveGradient',true, 'Display','iter','HessianApproximation','lbfgs')`.
- Lines 61-62: computes `target_fun` using `target_fun=@(x)tgrape(spin_system,drift,controls,waveform, x,time_unit,rho_init,rho_targ)`.

## Implementation structure

- Optimisation of slice durations in a composite inversion
- pulse with specified amplitudes, phases, and a constrain-
- ed overall duration. The initial guess is
- 270(-x)360(x)90(y)270(-y)360(y)90(x) [Fig. 3]
- from https://doi.org/10.1016/0022-2364(83)90133-6 --the
- optimisation demonstrates that a slightly better pulse of
- the same power and duration exists.
- Calculation time: minutes.
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals over the area
- that needs to be affected by the pulse (25 kHz either side)
- Select a basis set -IK-2 keeps complete basis on each

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimoptions()`, `tgrape()`, `fmincon()`, `step()`, `waveform()`, `dt_old()`, `dt_new()`, `liquid()`, `apodisation()`.
