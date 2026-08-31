# examples/optimal_control/state_transfer_wf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_wf.m`
- Signature: `state_transfer_wf()`
- Total lines: 92

## Purpose

A transfer of population from the lowermost energy level in a four-spin system to the uppermost energy level using the wave- function space version of the GRAPE algorithm. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 15-16: Isotopes; implemented by `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 18-19: Interactions; implemented by `inter.zeeman.scalar={1.5, 2.0, 30.0, 40.0}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='zeeman-wavef'`.
- Lines 30-31: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Bottom ground state to start; implemented by `rho_init=zeros(16,1); rho_init(end)=1`.
- Lines 37-38: Top excited state to finish; implemented by `rho_targ=zeros(16,1); rho_targ(1)=1`.
- Lines 40-41: Get the control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 46-47: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 49-50: Transmitter offsets; implemented by `parameters.spins={'1H','13C'}`.
- Lines 54-55: Define control parameters; implemented by `control.isotopes={'1H','13C'}`.
- Lines 68-69: Plots during optimisation; implemented by `control.plotting={'xy_controls','spectrogram'}`.
- Lines 71-72: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 74-75: Initial guess; implemented by `guess=randn(4,100)/10`.
- Lines 77-78: Run the optimisation, get normalised pulse; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 80-81: Apply power level scaling; implemented by `pulse=mean(control.pwr_levels)*pulse`.
- Lines 84-85: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5, 2.0, 30.0, 40.0}`.
- Lines 20: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 21: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 22: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=150`.
- Lines 23: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=150`.
- Lines 24: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=50`.
- Lines 27: computes `bas.formalism` using `bas.formalism='zeeman-wavef'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `rho_init` using `rho_init=zeros(16,1); rho_init(end)=1`.
- Lines 38: computes `rho_targ` using `rho_targ=zeros(16,1); rho_targ(1)=1`.
- Lines 41: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 42: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 43: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 44: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 47: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- A transfer of population from the lowermost energy level in a
- four-spin system to the uppermost energy level using the wave-
- function space version of the GRAPE algorithm.
- Calculation time: minutes.
- Magnetic field
- Isotopes
- Interactions
- Basis set
- Run Spinach housekeeping
- Bottom ground state to start
- Top excited state to finish
- Get the control operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rho_init()`, `rho_targ()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
