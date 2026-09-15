# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m`
- Signature: `mq_excitation()`
- Total lines: 108

## Purpose

Optimal control design of the multiple-quantum excitation pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the excitation pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is three rotor periods (240 us) long in 0.5 us slices, the controls are Cartesian, and the 100 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is Iz and the target is the Hermitian combination of the +MQ and -MQ coherences between the m=+3/2 and m=-3/2 levels (3Q) or between the m=+5/2 and m=-5/2 levels (5Q), as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation.

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

- Lines 27-28: Coherence order to excite, 3 or 5; implemented by `mq_order=5`.
- Lines 30-31: 400 MHz magnet; implemented by `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 34-35: Quadrupolar coupling and shielding anisotropy; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 39-40: Hilbert space formalism; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Rotor phase resolved drift Hamiltonians; implemented by `parameters.spins={'27Al'}`.
- Lines 57-58: Initial state, Iz; implemented by `rho_init=state(spin_system,'Lz','27Al')`.
- Lines 61-62: Target state, symmetric MQ coherence between m=+mq_order/2 and m=-mq_order/2; implemented by `rho_targ=zeros(6); rho_targ(3.5-mq_order/2,3.5+mq_order/2)=1`.
- Lines 65-66: Control operators; implemented by `Lx=operator(spin_system,'Lx','27Al')`.
- Lines 69-70: Control parameters; implemented by `control.isotopes={'27Al'}`.
- Lines 82-83: Plotting options; implemented by `control.plotting={'amp_controls','phi_controls','spectrogram'}`.
- Lines 85-86: Random initial guess, amplitudes up to 10% of the ceiling; implemented by `amp=0.1*rand(1,480); phi=2*pi*rand(1,480); amp(randi(480))=1`.
- Lines 89-90: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 92-93: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 95-96: Clip the amplitude to the ceiling; implemented by `[amp,phi]=cartesian2polar(pulse(1,:),pulse(2,:)); amp=min(amp,1)`.
- Lines 99-100: Report the fidelity of the clipped pulse; implemented by `[~,fidelity]=grape_xy(pulse,spin_system)`.
- Lines 103-104: Save the waveform in rad/s; implemented by `pulse=control.pwr_levels*pulse; pulse_dt=control.pulse_dt`.

### Key state/data transformations

- Lines 28: computes `mq_order` using `mq_order=5`.
- Lines 31: computes `sys.magnet` using `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 32: computes `sys.isotopes` using `sys.isotopes={'27Al'}`.
- Lines 35: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 36: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-5 -5 10]}`.
- Lines 37: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 40: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `parameters.spins` using `parameters.spins={'27Al'}`.
- Lines 50: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 51: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 52: computes `parameters.n_ticks` using `parameters.n_ticks=160`.
- Lines 53: computes `parameters.n_phases` using `parameters.n_phases=80`.
- Lines 54: computes `parameters.n_slices` using `parameters.n_slices=480`.
- Lines 55: computes `control.drifts` using `control.drifts=mqmas_drifts(spin_system,parameters)`.
- Lines 58: computes `rho_init` using `rho_init=state(spin_system,'Lz','27Al')`.
- Lines 62: computes `rho_targ` using `rho_targ=zeros(6); rho_targ(3.5-mq_order/2,3.5+mq_order/2)=1`.

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
