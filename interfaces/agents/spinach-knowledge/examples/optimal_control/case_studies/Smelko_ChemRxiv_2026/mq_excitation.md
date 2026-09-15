# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m`
- Signature: `mq_excitation()`
- Total lines: 110

## Purpose

Optimal control design of the multiple-quantum excitation pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi- nach, the excitation pulse optimisation from A single 27Al nucleus with the quadrupolar coupling and the shi- elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to secon

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Coherence order to excite, 3 or 5; implemented by `mq_order=5`.
- Lines 32-33: 400 MHz magnet; implemented by `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 36-37: Quadrupolar coupling and shielding anisotropy; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 41-42: Hilbert space formalism; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 45-46: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 50-51: Rotor phase resolved drift Hamiltonians; implemented by `parameters.spins={'27Al'}`.
- Lines 59-60: Initial state, Iz; implemented by `rho_init=state(spin_system,'Lz','27Al')`.
- Lines 63-64: Target state, symmetric MQ coherence between m=+mq_order/2 and m=-mq_order/2; implemented by `rho_targ=zeros(6); rho_targ(3.5-mq_order/2,3.5+mq_order/2)=1`.
- Lines 67-68: Control operators; implemented by `Lx=operator(spin_system,'Lx','27Al')`.
- Lines 71-72: Control parameters; implemented by `control.isotopes={'27Al'}`.
- Lines 84-85: Plotting options; implemented by `control.plotting={'amp_controls','phi_controls','spectrogram'}`.
- Lines 87-88: Random initial guess, amplitudes up to 10% of the ceiling; implemented by `amp=0.1*rand(1,480); phi=2*pi*rand(1,480); amp(randi(480))=1`.
- Lines 91-92: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 94-95: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 97-98: Clip the amplitude to the ceiling; implemented by `[amp,phi]=cartesian2polar(pulse(1,:),pulse(2,:)); amp=min(amp,1)`.
- Lines 101-102: Report the fidelity of the clipped pulse; implemented by `[~,fidelity]=grape_xy(pulse,spin_system)`.
- Lines 105-106: Save the waveform in rad/s; implemented by `pulse=control.pwr_levels*pulse; pulse_dt=control.pulse_dt`.

### Key state/data transformations

- Lines 30: computes `mq_order` using `mq_order=5`.
- Lines 33: computes `sys.magnet` using `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 34: computes `sys.isotopes` using `sys.isotopes={'27Al'}`.
- Lines 37: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 38: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-5 -5 10]}`.
- Lines 39: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 42: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 43: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 46: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 51: computes `parameters.spins` using `parameters.spins={'27Al'}`.
- Lines 52: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 53: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 54: computes `parameters.n_ticks` using `parameters.n_ticks=160`.
- Lines 55: computes `parameters.n_phases` using `parameters.n_phases=80`.
- Lines 56: computes `parameters.n_slices` using `parameters.n_slices=480`.
- Lines 57: computes `control.drifts` using `control.drifts=mqmas_drifts(spin_system,parameters)`.
- Lines 60: computes `rho_init` using `rho_init=state(spin_system,'Lz','27Al')`.
- Lines 64: computes `rho_targ` using `rho_targ=zeros(6); rho_targ(3.5-mq_order/2,3.5+mq_order/2)=1`.

## Implementation structure

- Optimal control design of the multiple-quantum excitation pulse
- of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
- nach, the excitation pulse optimisation from
- A single 27Al nucleus with the quadrupolar coupling and the shi-
- elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
- ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
- magnet. The quadrupolar interaction is taken to second order in
- the rotating frame, and the powder average runs over 200 crystal-
- lite orientations at 80 initial rotor phases each. The pulse is
- three rotor periods (240 us) long in 0.5 us slices, the controls
- are Cartesian, and the 100 kHz amplitude ceiling is enforced by
- a spillout penalty followed by clipping. The initial state is Iz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `assume()`, `mqmas_drifts()`, `state()`, `operator()`, `randi()`, `optimcon()`, `fmaxnewton()`, `cartesian2polar()`, `polar2cartesian()`, `grape_xy()`.
