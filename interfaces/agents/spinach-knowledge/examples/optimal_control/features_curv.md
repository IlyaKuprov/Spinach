# examples/optimal_control/features_curv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_curv.m`
- Signature: `features_curv()`
- Total lines: 105

## Purpose

A transfer of coherence from longitudinal magnetization into a two-spin singlet state with a distribution of B1 powers. An ensemble of ten spin systems with different power levels is simultaneously driven to optimal fidelity, which in this case is 1/sqrt(2) = 0.7071 Curvilinear GRAPE interface is used -the user specifies the definition of the curvilinear coordinates and the Jacobian. In this case, the coor- dinates a

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

- Lines 16-17: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Isotopes; implemented by `sys.isotopes={'13C','13C'}`.
- Lines 22-23: Interactions; implemented by `inter.zeeman.scalar={0.00 0.25}`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,'Lz','all')`.
- Lines 39-40: Set up and normalise the target state; implemented by `rho_targ=singlet(spin_system,1,2)`.
- Lines 43-44: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 47-48: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 50-51: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 64-66: Diagnostic output; implemented by `control.plotting={'correlation_order','coherence_order', 'xy_controls','robustness','spectrogram'}`.
- Lines 68-69: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 71-72: Curvilinear coordinates; implemented by `u2x=@(u)[u(1)*cos(u(2))`.
- Lines 75-76: Curvilinear coordinate Jacobian; implemented by `dx_du=@(u)[ cos(u(2)) sin(u(2))`.
- Lines 79-80: Initial guess; implemented by `guess_u=[ones(1,50); (pi/2)*randn(1,50)]`.
- Lines 82-83: Run LBFGS GRAPE pulse optimisation; implemented by `curv_profile=fmaxnewton(spin_system,@(x,y)grape_curv(x,u2x,dx_du,y),guess_u)`.
- Lines 85-86: Get Cartesian components of the pulse; implemented by `cart_profile=zeros(size(curv_profile))`.
- Lines 93-95: Simulate the optimised pulse; implemented by `rho=shaped_pulse_xy(spin_system,D,{Lx,Ly},{CLx,CLy}, control.pulse_dt,rho_init,'expv-pwc')`.

### Control flow inferred from the code

- Line 87: `for` loop over `n=1:size(curv_profile,2)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'13C','13C'}`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.25}`.
- Lines 24: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 25: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=60.0`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `rho_init` using `rho_init=state(spin_system,'Lz','all')`.
- Lines 40: computes `rho_targ` using `rho_targ=singlet(spin_system,1,2)`.
- Lines 44: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 45: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 48: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 51: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 52: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 53: computes `control.drifts` using `control.drifts={{D}}`.
- Lines 54: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 55: computes `control.rho_init` using `control.rho_init={rho_init}`.

## Implementation structure

- A transfer of coherence from longitudinal magnetization into a two-spin
- singlet state with a distribution of B1 powers. An ensemble of ten spin
- systems with different power levels is simultaneously driven to optimal
- fidelity, which in this case is 1/sqrt(2) = 0.7071
- Curvilinear GRAPE interface is used -the user specifies the definition
- of the curvilinear coordinates and the Jacobian. In this case, the coor-
- dinates are phase-amplitude.
- Calculation time: minutes.
- Magnetic field
- Isotopes
- Interactions
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `singlet()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `grape_curv()`, `cart_profile()`, `u2x()`, `curv_profile()`, `shaped_pulse_xy()`, `coherence()`, `correlation()`.
