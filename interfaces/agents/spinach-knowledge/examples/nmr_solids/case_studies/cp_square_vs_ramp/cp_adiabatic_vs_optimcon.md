# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m`
- Signature: `cp_adiabatic_vs_optimcon()`
- Total lines: 147

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame using (a) tangent-ramped adiabatic CP; (b) numerically optimised (GRAPE method) shortcut to adiabaticity. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=9.394`.
- Lines 15-16: Interactions; implemented by `inter.zeeman.scalar={0.00 0.00}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-32: % Tangent ramp CP simulation; implemented by `parameters.time_steps=2e-6*ones(1,500)`.
- Lines 31-32: Common experiment parameters; implemented by `parameters.time_steps=2e-6*ones(1,500)`.
- Lines 41-42: Simulate tangent ramped amplitude CP; implemented by `ramp_up=tan(linspace(-1.4,1.4,500))`.
- Lines 50-51: Plotting -waveform; implemented by `kfigure(); scale_figure([1.5 1.5])`.
- Lines 59-60: Plotting -trajectory; implemented by `subplot(2,2,4); plot(time_axis,real(fid_a)); hold on`.
- Lines 64-67: % GRAPE optimisation with the same timing; implemented by `control.isotopes={'1H','15N'}`.
- Lines 66-67: Drift Hamiltonians for every system in the powder; implemented by `control.isotopes={'1H','15N'}`.
- Lines 71-72: Initial state: Ly on 1H; implemented by `rho_init=state(spin_system,'Ly','1H')`.
- Lines 76-77: Target state: Lx on 15N; implemented by `rho_targ=state(spin_system,'Lx','15N')`.
- Lines 81-82: Control operators; implemented by `LyH=operator(spin_system,'Ly','1H')`.
- Lines 86-87: Other GRAPE settings; implemented by `control.pwr_levels=2*pi*5e4`.
- Lines 94-95: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 97-98: Initial guess from above; implemented by `guess=[ramp_down; ramp_up]`.
- Lines 100-101: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.00}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 19: computes `inter.temperature` using `inter.temperature=298`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.time_steps` using `parameters.time_steps=2e-6*ones(1,500)`.
- Lines 33-34: computes `parameters.irr_opers` using `parameters.irr_opers={operator(spin_system,'Ly','1H') operator(spin_system,'Lx','15N')}`.
- Lines 35: computes `parameters.exc_opers` using `parameters.exc_opers={operator(spin_system,'Lx','1H')}`.
- Lines 36: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','15N')`.
- Lines 37: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 38: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'15N'}`.
- Lines 42: computes `ramp_up` using `ramp_up=tan(linspace(-1.4,1.4,500))`.
- Lines 45: computes `ramp_down` using `ramp_down=fliplr(ramp_up)`.
- Lines 46: computes `irr_powers_a` using `irr_powers_a=[5e4*ramp_down; 5e4*ramp_up]`.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) tangent-ramped adiabatic CP; (b) numerically
- optimised (GRAPE method) shortcut to adiabaticity.
- Calculation time: minutes
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- % Tangent ramp CP simulation
- Common experiment parameters
- Simulate tangent ramped amplitude CP
- Plotting -waveform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `fliplr()`, `powder()`, `kfigure()`, `scale_figure()`, `cumsum()`, `subplot()`, `time_axis()`, `kxlabel()`, `kylabel()`, `ylim()`, `klegend()`, `drifts()`.
