# examples/optimal_control/pattern_pulse_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/pattern_pulse_2.m`
- Signature: `pattern_pulse_2()`
- Total lines: 106

## Purpose

Transmitter offset selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified transmitter offset intervals have magnetisation arriving into us- er specified states. The pulse is phase-modulated. Calculation time: minutes.

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
- Lines 18-19: Put the spin at 0 ppm; implemented by `inter.zeeman.scalar{1}=0`.
- Lines 21-22: No approximations; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Get pertinent spin states; implemented by `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(Sx,2)`.
- Lines 33-34: Get pertinent control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 38-39: Get the drift Hamiltonian (zero here); implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 41-42: Transmitter offset range in Hz; implemented by `toff_range=linspace(-4e3,4e3,128)`.
- Lines 44-45: Excitation pattern; implemented by `X=ones(1,128); X(1:20)=0; X(109:128)=0; X(54:73)=0`.
- Lines 52-53: Patterened target state array; implemented by `rho_targ=X.*repmat(Sx,1,128)+Z.*repmat(Sz,1,128)`.
- Lines 55-56: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 71-73: Plotting options; implemented by `control.plotting={'phi_controls','xy_controls', 'robustness','spectrogram'}`.
- Lines 75-76: Initial guess; implemented by `guess=(pi/4)*ones(1,300)`.
- Lines 78-79: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 81-82: Run LBFGS GRAPE pulse optimisation; implemented by `phi_profile=fmaxnewton(spin_system,@grape_phase,guess)`.
- Lines 84-85: Loop over the offsets; implemented by `X=zeros(1,128); Z=zeros(1,128)`.
- Lines 88-89: Get Cartesian components of the pulse; implemented by `amp_profile=control.pwr_levels*control.amplitudes`.

### Control flow inferred from the code

- Line 86: `parfor` loop over `n=1:numel(toff_range)`.

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
- Lines 36: computes `Lz` using `Lz=operator(spin_system,'Lz','13C')`.
- Lines 39: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 42: computes `toff_range` using `toff_range=linspace(-4e3,4e3,128)`.
- Lines 45: computes `X` using `X=ones(1,128); X(1:20)=0; X(109:128)=0; X(54:73)=0`.
- Lines 46: computes `Z` using `Z=zeros(1,128); Z(1:20)=1; Z(109:128)=1; Z(54:73)=1`.
- Lines 53: computes `rho_targ` using `rho_targ=X.*repmat(Sx,1,128)+Z.*repmat(Sz,1,128)`.
- Lines 56: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 57: computes `control.channels` using `control.channels=[1; 1]`.

## Implementation structure

- Transmitter offset selective excitation described in Glaser group
- paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified
- transmitter offset intervals have magnetisation arriving into us-
- er specified states. The pulse is phase-modulated.
- Calculation time: minutes.
- Magnetic field
- Single carbon spin
- Put the spin at 0 ppm
- No approximations
- Run Spinach housekeeping
- Get pertinent spin states
- Get pertinent control operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `figure()`, `kxlabel()`, `klegend()`, `ylim()`, `num2cell()`, `optimcon()`, `fmaxnewton()`, `polar2cartesian()`, `shaped_pulse_xy()`, `toff_range()`.
