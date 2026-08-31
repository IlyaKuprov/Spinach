# examples/optimal_control/state_transfer_m2s.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_m2s.m`
- Signature: `state_transfer_m2s()`
- Total lines: 147

## Purpose

A transfer of coherence from longitudinal magnetization into a two-spin singlet state in allyl pyruvate with a distribution of B1 powers and transmitter offsets. XX and YY components of the singlet dephase rapidly in this system, and are therefore drop- ped from the target state specification. Calculation time: many hours.

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

- Lines 14-15: Get the spin system; implemented by `[sys,inter]=allyl_pyruvate({'1H'})`.
- Lines 17-18: Kill the methyl group; implemented by `sys.isotopes=sys.isotopes(1:5)`.
- Lines 25-26: Magnetic field (500.13 MHz); implemented by `sys.magnet=11.7464`.
- Lines 28-29: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,'Lz','all')`.
- Lines 39-42: Set up and normalise target state; implemented by `rho_targ=state(spin_system,{'Lz','Lz'}, {idxof(sys,'Hb'), idxof(sys,'Hc')})`.
- Lines 44-45: Get the control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 48-49: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 51-52: Transmitter offset; implemented by `parameters.spins={'1H'}; parameters.offset=2670`.
- Lines 55-56: Offset distribution generator; implemented by `Hz=operator(spin_system,'Lz','1H')`.
- Lines 58-59: Define control parameters; implemented by `control.isotopes={'1H'}`.
- Lines 74-77: Plots during optimisation; implemented by `control.plotting={'correlation_order','coherence_order', 'xy_controls','local_each_spin', 'amp_controls','spectrogram'}`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 82-83: Initial guess -smack all spins; implemented by `time_axis=cumsum(control.pulse_dt)`.
- Lines 87-88: Run the optimisation, get normalised pulse; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,pulse)`.
- Lines 90-93: % Simple pulse-acquire; implemented by `parameters.rho0=rho_init`.
- Lines 92-93: Parameters; implemented by `parameters.rho0=rho_init`.

### Key state/data transformations

- Lines 15: computes `[sys,inter]` using `[sys,inter]=allyl_pyruvate({'1H'})`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes=sys.isotopes(1:5)`.
- Lines 19: computes `sys.labels` using `sys.labels=sys.labels(1:5)`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=inter.zeeman.scalar(1:5)`.
- Lines 21: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=inter.zeeman.matrix(1:5)`.
- Lines 22: computes `inter.coordinates` using `inter.coordinates=inter.coordinates(1:5)`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=inter.coupling.scalar(1:5,1:5)`.
- Lines 26: computes `sys.magnet` using `sys.magnet=11.7464`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `rho_init` using `rho_init=state(spin_system,'Lz','all')`.
- Lines 40-42: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz','Lz'}, {idxof(sys,'Hb'), idxof(sys,'Hc')})`.
- Lines 45: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 46: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 49: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 52: computes `parameters.spins` using `parameters.spins={'1H'}; parameters.offset=2670`.
- Lines 56: computes `Hz` using `Hz=operator(spin_system,'Lz','1H')`.

## Implementation structure

- A transfer of coherence from longitudinal magnetization into a
- two-spin singlet state in allyl pyruvate with a distribution of
- B1 powers and transmitter offsets. XX and YY components of the
- singlet dephase rapidly in this system, and are therefore drop-
- ped from the target state specification.
- Calculation time: many hours.
- Get the spin system
- Kill the methyl group
- Magnetic field (500.13 MHz)
- Formalism and basis set
- Run Spinach housekeeping
- Set up and normalise the initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `allyl_pyruvate()`, `create()`, `basis()`, `state()`, `idxof()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `cumsum()`, `fmaxnewton()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`.
