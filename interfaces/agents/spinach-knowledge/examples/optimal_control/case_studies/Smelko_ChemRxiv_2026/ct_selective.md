# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m`
- Signature: `ct_selective()`
- Total lines: 102

## Purpose

Optimal control design of the central transition selective pulse of the z-filtered 27Al MQMAS experiment. Reproduces, using Spinach, the soft pulse optimisation from https://doi.org/10.26434/chemrxiv.15008427 A single 27Al nucleus with the quadrupolar coupling and the shielding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10 ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz magnet. The quadrupolar interaction is taken to second order in the rotating frame, and the powder average runs over 200 crystallite orientations at 80 initial rotor phases each. The pulse is 50 us long in 0.5 us slices, the controls are Cartesian, and the 10 kHz amplitude ceiling is enforced by a spillout penalty followed by clipping. The initial state is the population difference across the central transition, and the target is the single-quantum coherence of the central transition, as in the paper. The resulting waveform is saved for the MQMAS efficiency calculation.

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

- Lines 26-27: 400 MHz magnet; implemented by `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 30-31: Quadrupolar coupling and shielding anisotropy; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 35-36: Hilbert space formalism; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Rotor phase resolved drift Hamiltonians; implemented by `parameters.spins={'27Al'}`.
- Lines 53-54: Initial state, population difference across the central transition; implemented by `rho_init=diag([0 0 1 -1 0 0]); rho_init=rho_init/norm(rho_init,'fro')`.
- Lines 56-57: Target state, single-quantum coherence of the central transition; implemented by `rho_targ=zeros(6); rho_targ(3,4)=1`.
- Lines 59-60: Control operators; implemented by `Lx=operator(spin_system,'Lx','27Al')`.
- Lines 63-64: Control parameters; implemented by `control.isotopes={'27Al'}`.
- Lines 76-77: Plotting options; implemented by `control.plotting={'amp_controls','phi_controls','spectrogram'}`.
- Lines 79-80: Random initial guess, amplitudes up to 10% of the ceiling with one slice at the ceiling; implemented by `amp=0.1*rand(1,100); phi=2*pi*rand(1,100); amp(randi(100))=1`.
- Lines 83-84: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 86-87: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 89-90: Clip the amplitude to the ceiling; implemented by `[amp,phi]=cartesian2polar(pulse(1,:),pulse(2,:)); amp=min(amp,1)`.
- Lines 93-94: Report the fidelity of the clipped pulse; implemented by `[~,fidelity]=grape_xy(pulse,spin_system)`.
- Lines 97-98: Save the waveform in rad/s; implemented by `pulse=control.pwr_levels*pulse; pulse_dt=control.pulse_dt`.

### Key state/data transformations

- Lines 27: computes `sys.magnet` using `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'27Al'}`.
- Lines 31: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0])`.
- Lines 32: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-5 -5 10]}`.
- Lines 33: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 36: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'27Al'}`.
- Lines 46: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 47: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 48: computes `parameters.n_ticks` using `parameters.n_ticks=160`.
- Lines 49: computes `parameters.n_phases` using `parameters.n_phases=80`.
- Lines 50: computes `parameters.n_slices` using `parameters.n_slices=100`.
- Lines 51: computes `control.drifts` using `control.drifts=mqmas_drifts(spin_system,parameters)`.
- Lines 54: computes `rho_init` using `rho_init=diag([0 0 1 -1 0 0]); rho_init=rho_init/norm(rho_init,'fro')`.
- Lines 57: computes `rho_targ` using `rho_targ=zeros(6); rho_targ(3,4)=1`.
- Lines 60: computes `Lx` using `Lx=operator(spin_system,'Lx','27Al')`.

## Implementation structure

- Optimal control design of the central transition selective pulse
- of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
- nach, the soft pulse optimisation from
- A single 27Al nucleus with the quadrupolar coupling and the shi-
- elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
- ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
- magnet. The quadrupolar interaction is taken to second order in
- the rotating frame, and the powder average runs over 200 crystal-
- lite orientations at 80 initial rotor phases each. The pulse is
- 50 us long in 0.5 us slices, the controls are Cartesian, and the
- 10 kHz amplitude ceiling is enforced by a spillout penalty follo-
- wed by clipping. The initial state is the population difference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `assume()`, `mqmas_drifts()`, `operator()`, `randi()`, `optimcon()`, `fmaxnewton()`, `cartesian2polar()`, `polar2cartesian()`, `grape_xy()`.
