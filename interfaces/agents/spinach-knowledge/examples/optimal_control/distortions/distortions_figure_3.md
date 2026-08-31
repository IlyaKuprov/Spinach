# examples/optimal_control/distortions/distortions_figure_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_3.m`
- Signature: `distortions_figure_3()`
- Total lines: 225

## Purpose

Figure 3 from the paper by Rasulov and Kuprov:

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

- Lines 10-11: Set the magnetic field; implemented by `sys.magnet=28.18`.
- Lines 13-15: Put 100 non-interacting spins at equal intervals within the [-100,+100] ppm chemical shift range; implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 21-23: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Set up spin states; implemented by `Sx=state(spin_system,'Lx','13C')`.
- Lines 40-41: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 44-45: Get the drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 47-48: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 61-62: Last five slices are dead time; implemented by `control.freeze=zeros(2,45)`.
- Lines 65-66: Plotting options; implemented by `control.plotting={'xy_controls','robustness','spectrogram'}`.
- Lines 68-69: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 71-72: Initial guess; implemented by `guess=(1/4)*ones(2,45)`.
- Lines 75-76: GRAPE optimisation; implemented by `xy_profile=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 78-79: Extract the waveform; implemented by `xy_profile=mean(control.pwr_levels)*xy_profile`.
- Lines 82-83: Simulate the optimised pulse; implemented by `rho_init=state(spin_system,'Lz','13C')`.
- Lines 87-88: Set acquisition parameters; implemented by `parameters.spins={'13C'}`.
- Lines 99-100: Simulate the free induction decay; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 102-103: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.

### Control flow inferred from the code

- Line 16: `for` loop over `n=1:n_spins`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=28.18`.
- Lines 15: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 17: computes `sys.isotopes{n}` using `sys.isotopes{n}='13C'`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins))`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 25: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `Sx` using `Sx=state(spin_system,'Lx','13C')`.
- Lines 34: computes `Sy` using `Sy=state(spin_system,'Ly','13C')`.
- Lines 35: computes `Sz` using `Sz=state(spin_system,'Lz','13C')`.
- Lines 41: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 42: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 45: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 48: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 49: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 50: computes `control.drifts` using `control.drifts={{D}}`.

## Implementation structure

- Figure 3 from the paper by Rasulov and Kuprov:
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping
- Set up spin states
- Get the control operators
- Get the drift Hamiltonian
- Define control parameters
- Last five slices are dead time

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `guess()`, `fmaxnewton()`, `xy_profile()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `figure()`.
