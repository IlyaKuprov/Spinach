# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m`
- Signature: `bb_inversion_pulse()`
- Total lines: 101

## Purpose

Broadband inversion pulse design for liquid-state NMR. Reprodu- ces, using Spinach, the second example from: A single proton is considered in the rotating frame with multiple transmitter offsets (or chemical shifts). The goal is to design a 600 µs broadband inversion pulse (1 µs slices) that performs: I_z → -I_z uniformly over a frequency offset range of ±50 kHz; controls are Cartesian (Lx, Ly) operators in the rotat

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

- Lines 19-20: Magnetic field (Tesla); implemented by `sys.magnet=14.1`.
- Lines 23-24: Chemical shift (ppm); implemented by `inter.zeeman.scalar={0.0}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Initial and target states; implemented by `Sz=state(spin_system,'Lz',1)`.
- Lines 38-39: Control and offset operators; implemented by `LxH=operator(spin_system,'Lx',1)`.
- Lines 43-44: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 46-47: Control data structure; implemented by `control.isotopes={'1H'}`.
- Lines 64-65: Random guess; implemented by `guess=randn(2,600)/10`.
- Lines 67-68: Optimisation; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 71-72: Return to physical units; implemented by `rf_scale=mean(control.pwr_levels)`.
- Lines 76-77: Offset grid for verification; implemented by `offs_hz=linspace(-100e3,100e3,201)`.
- Lines 80-81: Test simulation; implemented by `parfor k=1:numel(offs_hz)`.
- Lines 83-84: Add the offset term; implemented by `Hd=H+2*pi*offs_hz(k)*LzH`.
- Lines 86-88: Run the pulse; implemented by `rho_f=shaped_pulse_xy(spin_system,Hd,{LxH, LyH},{CLx,CLy}, control.pulse_dt,Sz,'expv-pwc')`.
- Lines 90-91: Compute inversion efficiency; implemented by `inv_eff(k)=-real(Sz'*rho_f)`.
- Lines 95-96: Plot the fidelity profile; implemented by `kfigure(); plot(offs_hz/1e3,inv_eff); kgrid`.

### Control flow inferred from the code

- Line 81: `parfor` loop over `k=1:numel(offs_hz)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `Sz` using `Sz=state(spin_system,'Lz',1)`.
- Lines 39: computes `LxH` using `LxH=operator(spin_system,'Lx',1)`.
- Lines 40: computes `LyH` using `LyH=operator(spin_system,'Ly',1)`.
- Lines 41: computes `LzH` using `LzH=operator(spin_system,'Lz',1)`.
- Lines 44: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 47: computes `control.isotopes` using `control.isotopes={'1H'}`.
- Lines 48: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 49: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 50: computes `control.operators` using `control.operators={LxH,LyH}`.
- Lines 51: computes `control.off_ops` using `control.off_ops={LzH}`.
- Lines 52: computes `control.offsets` using `control.offsets={linspace(-50e3,50e3,101)}`.
- Lines 53: computes `control.rho_init` using `control.rho_init={+Sz}`.

## Implementation structure

- Broadband inversion pulse design for liquid-state NMR. Reprodu-
- ces, using Spinach, the second example from:
- A single proton is considered in the rotating frame with multiple
- transmitter offsets (or chemical shifts). The goal is to design a
- 600 µs broadband inversion pulse (1 µs slices) that performs:
- I_z → -I_z
- uniformly over a frequency offset range of ±50 kHz; controls are
- Cartesian (Lx, Ly) operators in the rotating frame.
- Magnetic field (Tesla)
- Chemical shift (ppm)
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `xy_profile()`, `offs_hz()`, `shaped_pulse_xy()`, `inv_eff()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`.
