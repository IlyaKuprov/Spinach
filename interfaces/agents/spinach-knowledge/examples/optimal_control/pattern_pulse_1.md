# examples/optimal_control/pattern_pulse_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/pattern_pulse_1.m`
- Signature: `pattern_pulse_1()`
- Total lines: 103

## Purpose

Nutation frequency selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified nutation frequency intervals have magnetisation arriving into us- er specified states. The pulse is phase-modulated. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=28.18`.
- Lines 15-16: Single carbon spin; implemented by `sys.isotopes{1}='13C'`.
- Lines 18-19: Transmitter is on resonance; implemented by `inter.zeeman.scalar{1}=0`.
- Lines 21-22: No approximations; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Get pertinent spin states; implemented by `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(Sx,2)`.
- Lines 33-34: Get pertinent control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 37-38: Get the drift Hamiltonian (zero here); implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 40-41: Nutation frequency range in Hz; implemented by `nutf_range=linspace(6e3,14e3,128)`.
- Lines 43-44: Excitation pattern; implemented by `X=ones(1,128); X(1:20)=0; X(109:128)=0; X(54:73)=0`.
- Lines 51-52: Patterened target state array; implemented by `rho_targ=X.*repmat(Sx,1,128)+Z.*repmat(Sz,1,128)`.
- Lines 54-55: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 68-70: Plotting options; implemented by `control.plotting={'phi_controls','xy_controls', 'robustness','spectrogram'}`.
- Lines 72-73: Initial guess; implemented by `guess=(pi/4)*ones(1,250)`.
- Lines 75-76: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 78-79: Run LBFGS GRAPE pulse optimisation; implemented by `phi_profile=fmaxnewton(spin_system,@grape_phase,guess)`.
- Lines 81-82: Loop over the power levels; implemented by `X=zeros(1,128); Z=zeros(1,128)`.
- Lines 85-86: Get Cartesian components of the pulse; implemented by `amp_profile=control.pwr_levels(n)*control.amplitudes`.

### Control flow inferred from the code

- Line 83: `parfor` loop over `n=1:numel(control.pwr_levels)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=28.18`.
- Lines 16: computes `sys.isotopes{1}` using `sys.isotopes{1}='13C'`.
- Lines 19: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `Sx` using `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(Sx,2)`.
- Lines 31: computes `Sz` using `Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(Sz,2)`.
- Lines 34: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 35: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 38: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 41: computes `nutf_range` using `nutf_range=linspace(6e3,14e3,128)`.
- Lines 44: computes `X` using `X=ones(1,128); X(1:20)=0; X(109:128)=0; X(54:73)=0`.
- Lines 45: computes `Z` using `Z=zeros(1,128); Z(1:20)=1; Z(109:128)=1; Z(54:73)=1`.
- Lines 52: computes `rho_targ` using `rho_targ=X.*repmat(Sx,1,128)+Z.*repmat(Sz,1,128)`.
- Lines 55: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 56: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 57: computes `control.drifts` using `control.drifts={{D}}`.

## Implementation structure

- Nutation frequency selective excitation described in Glaser group
- paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified
- nutation frequency intervals have magnetisation arriving into us-
- er specified states. The pulse is phase-modulated.
- Calculation time: minutes.
- Magnetic field
- Single carbon spin
- Transmitter is on resonance
- No approximations
- Run Spinach housekeeping
- Get pertinent spin states
- Get pertinent control operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `figure()`, `kxlabel()`, `klegend()`, `ylim()`, `num2cell()`, `optimcon()`, `fmaxnewton()`, `polar2cartesian()`, `shaped_pulse_xy()`.
