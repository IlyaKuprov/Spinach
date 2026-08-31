# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m`
- Signature: `bb_refocusing_pulse()`
- Total lines: 97

## Purpose

Spinach implementation of the broadband refocusing example from GRAPE is used to design a 200 µs broadband x-phase π pulse: {Sx -> Sx, Sy -> -Sy, Sz -> -Sz} over an offset range of ±12.5 kHz.

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

- Lines 15-16: Magnetic field (Tesla); implemented by `sys.magnet=14.1`.
- Lines 19-20: Chemical shift (ppm); implemented by `inter.zeeman.scalar={0.0}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Normalised Cartesian basis states; implemented by `Sx=state(spin_system,'Lx',1); Sx=Sx/norm(full(Sx),2)`.
- Lines 35-36: RF controls and offset operator; implemented by `Lx=operator(spin_system,'Lx',1)`.
- Lines 40-41: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 43-44: Control data structure; implemented by `control.isotopes={'1H'}`.
- Lines 59-61: Visual diagnostics; implemented by `control.plotting={'phi_controls','xy_controls', 'spectrogram','robustness'}`.
- Lines 63-64: Random initial guess; implemented by `guess=randn(2,600)/10`.
- Lines 66-67: Optimisation; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 70-71: Convert normalised waveform to physical rad/s controls; implemented by `CLx=control.pwr_levels*xy_profile(1,:)`.
- Lines 74-75: Test the pulse; implemented by `offs_hz=linspace(-25e3,25e3,201)`.
- Lines 79-80: Get drift Hamiltonian; implemented by `Hk=H+2*pi*offs_hz(k)*Lz`.
- Lines 82-84: Apply the pulse; implemented by `rho_k=shaped_pulse_xy(spin_system,Hk,control.operators,{CLx,CLy}, control.pulse_dt,[Sx,Sy,Sz],'expv-pwc')`.
- Lines 86-87: Calculate the fidelity; implemented by `fidelities(k)=real(trace([Sx,-Sy,-Sz]'*rho_k))/3`.
- Lines 91-92: Plot the fidelity profile; implemented by `kfigure(); plot(offs_hz/1e3,fidelities); kgrid`.

### Control flow inferred from the code

- Line 77: `parfor` loop over `k=1:numel(offs_hz)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `Sx` using `Sx=state(spin_system,'Lx',1); Sx=Sx/norm(full(Sx),2)`.
- Lines 32: computes `Sy` using `Sy=state(spin_system,'Ly',1); Sy=Sy/norm(full(Sy),2)`.
- Lines 33: computes `Sz` using `Sz=state(spin_system,'Lz',1); Sz=Sz/norm(full(Sz),2)`.
- Lines 36: computes `Lx` using `Lx=operator(spin_system,'Lx',1)`.
- Lines 37: computes `Ly` using `Ly=operator(spin_system,'Ly',1)`.
- Lines 38: computes `Lz` using `Lz=operator(spin_system,'Lz',1)`.
- Lines 41: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 44: computes `control.isotopes` using `control.isotopes={'1H'}`.
- Lines 45: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 46: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 47: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 48: computes `control.off_ops` using `control.off_ops={Lz}`.

## Implementation structure

- Spinach implementation of the broadband refocusing example from
- GRAPE is used to design a 200 µs broadband x-phase π pulse:
- {Sx -> Sx, Sy -> -Sy, Sz -> -Sz}
- over an offset range of ±12.5 kHz.
- Magnetic field (Tesla)
- Chemical shift (ppm)
- Basis set
- Spinach housekeeping
- Normalised Cartesian basis states
- RF controls and offset operator
- Drift Hamiltonian
- Control data structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `xy_profile()`, `offs_hz()`, `shaped_pulse_xy()`, `fidelities()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`.
