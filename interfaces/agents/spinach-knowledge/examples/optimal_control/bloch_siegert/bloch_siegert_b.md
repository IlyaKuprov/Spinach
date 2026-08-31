# examples/optimal_control/bloch_siegert/bloch_siegert_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/bloch_siegert_b.m`
- Signature: `bloch_siegert_b()`
- Total lines: 108

## Purpose

Bloch-Siegert shift compensation functionality demo. The script optimises a universal rotation pulse for a range of resonance offsets. As the control power is increased, Bloch-Siegert shift starts to reduce the fidelity unless it is correctly accounted for. Calculation time: minutes.

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

- Lines 13-14: Magnet field; implemented by `sys.magnet=28.18`.
- Lines 16-18: 100 non-interacting spins at equal intervals within [-100,+100] ppm chemical shift range; implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 24-26: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Set up spin states; implemented by `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 40-41: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 44-45: Get the drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 47-48: Optimal control settings; implemented by `control.isotopes={'13C'}`.
- Lines 59-60: Power levels to sweep (rad/s); implemented by `zeeman_frq=abs(spin('13C')*sys.magnet)`.
- Lines 63-64: Preallocate results; implemented by `fid_a=zeros(size(pwr_list))`.
- Lines 67-68: Common initial guess; implemented by `guess=randn(2,50)/10`.
- Lines 70-71: Loop over control power; implemented by `for n=1:numel(pwr_list)`.
- Lines 73-74: Set current power; implemented by `control.pwr_levels=pwr_list(n)`.
- Lines 76-77: Set pulse slice duration; implemented by `dt=(pi/100)/control.pwr_levels`.
- Lines 80-81: BSS settings: off, on; implemented by `control.bsiegert=false()`.
- Lines 86-87: Optimisation with and without BSS; implemented by `pulse_a=fmaxnewton(setting_a,@grape_xy,guess)`.
- Lines 90-91: Evaluation for a system with BSS; implemented by `[~,fid_a(n)]=ensemble(pulse_a,setting_b)`.
- Lines 96-97: Compute relative control powers; implemented by `relative_power=pwr_list./zeeman_frq`.

### Control flow inferred from the code

- Line 19: `for` loop over `n=1:n_spins`.
- Line 71: `for` loop over `n=1:numel(pwr_list)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=28.18`.
- Lines 18: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 20: computes `sys.isotopes{n}` using `sys.isotopes{n}='13C'`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins))`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 28: computes `bas.space_level` using `bas.space_level=1`.
- Lines 29: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `Sx` using `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 37: computes `Sy` using `Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2)`.
- Lines 38: computes `Sz` using `Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2)`.
- Lines 41: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 42: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 45: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 48: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 49: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 50: computes `control.drifts` using `control.drifts={{D}}`.

## Implementation structure

- Bloch-Siegert shift compensation functionality demo. The
- script optimises a universal rotation pulse for a range
- of resonance offsets. As the control power is increased,
- Bloch-Siegert shift starts to reduce the fidelity unless
- it is correctly accounted for.
- Calculation time: minutes.
- Magnet field
- 100 non-interacting spins at equal intervals
- within [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `spin()`, `pwr_list()`, `false()`, `optimcon()`, `true()`, `fmaxnewton()`, `fid_a()`, `ensemble()`, `fid_b()`.
