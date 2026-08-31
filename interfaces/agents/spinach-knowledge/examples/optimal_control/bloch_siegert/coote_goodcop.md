# examples/optimal_control/bloch_siegert/coote_goodcop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/coote_goodcop.m`
- Signature: `coote_goodcop()`
- Total lines: 147

## Purpose

Reproduction of the GOODCOP pulse design logic from Coote et al. with Bloch-Siegert corrections enabled in the optimiser and simulator The pulse enforces contracted-time C-alpha evolution while inverting CO

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `ppm2hz()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnetic field corresponding to 800 MHz 1H; implemented by `sys.magnet=18.8`.
- Lines 12-13: Single-spin carbon model; implemented by `sys.isotopes={'13C'}`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 20-21: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 24-25: Relevant operators and states; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 32-33: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 35-36: Paper parameters for GOODCOP; implemented by `carrier_ppm=100`.
- Lines 41-42: Offset grids from the paper; implemented by `ca_ppm=linspace(35,75,100)`.
- Lines 45-46: Convert to offset frequencies; implemented by `ca_hz=ppm2hz(sys.magnet,'13C',ca_ppm,carrier_ppm)`.
- Lines 50-51: Build state pairs correlated to offset ensemble; implemented by `rho_init=cell(1,numel(all_hz))`.
- Lines 67-68: Set up optimal control problem; implemented by `control.isotopes={'13C'}`.
- Lines 83-84: Enable Bloch-Siegert corrections; implemented by `control.bsiegert=true()`.
- Lines 86-87: Spinach housekeeping for Bloch-Siegert case; implemented by `spin_system_bs=optimcon(spin_system,control)`.
- Lines 89-90: Optimise the Bloch-Siegert waveform; implemented by `guess=randn(2,numel(control.pulse_dt))/8`.
- Lines 95-96: Disable Bloch-Siegert corrections; implemented by `control.bsiegert=false()`.
- Lines 98-99: Spinach housekeeping for no Bloch-Siegert case; implemented by `spin_system_nobs=optimcon(spin_system,control)`.
- Lines 101-102: Optimise the no Bloch-Siegert waveform; implemented by `pulse_nobs=fmaxnewton(spin_system_nobs,@grape_xy,guess)`.
- Lines 106-107: Validate on a dense offset grid; implemented by `eval_ppm=linspace(0,200,251)`.

### Control flow inferred from the code

- Line 53: `for` loop over `n=1:numel(ca_hz)`.
- Line 54: conditional branch on `mod(n,2)==1`.
- Line 61: `for` loop over `n=1:numel(co_hz)`.
- Line 111: `for` loop over `n=1:numel(eval_hz)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=18.8`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'13C'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 25: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 26: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 27: computes `Lz` using `Lz=operator(spin_system,'Lz','13C')`.
- Lines 28: computes `Ix` using `Ix=state(spin_system,'Lx','13C'); Ix=Ix/norm(full(Ix),2)`.
- Lines 29: computes `Iy` using `Iy=state(spin_system,'Ly','13C'); Iy=Iy/norm(full(Iy),2)`.
- Lines 30: computes `Iz` using `Iz=state(spin_system,'Lz','13C'); Iz=Iz/norm(full(Iz),2)`.
- Lines 33: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 36: computes `carrier_ppm` using `carrier_ppm=100`.
- Lines 37: computes `alpha_scale` using `alpha_scale=0.90`.
- Lines 38: computes `pulse_dur` using `pulse_dur=150e-6`.
- Lines 39: computes `max_rf_hz` using `max_rf_hz=15e3`.
- Lines 42: computes `ca_ppm` using `ca_ppm=linspace(35,75,100)`.

### Local helper functions

- Line 139: `ppm2hz()` — `function off_hz=ppm2hz(magnet,isotope,ppm_grid,carrier_ppm)`. Convert ppm values into offset frequencies in Hz
  - Representative operation: `frq_hz=abs(spin(isotope)*magnet/(2*pi))`.
  - Representative operation: `off_hz=(ppm_grid-carrier_ppm)*1e-6*frq_hz`.

## Implementation structure

- Reproduction of the GOODCOP pulse design logic from Coote et al.
- with Bloch-Siegert corrections enabled in the optimiser and simulator
- The pulse enforces contracted-time C-alpha evolution while inverting CO
- Magnetic field corresponding to 800 MHz 1H
- Single-spin carbon model
- Basis set
- Spinach housekeeping
- Relevant operators and states
- Drift Hamiltonian
- Paper parameters for GOODCOP
- Offset grids from the paper
- Convert to offset frequencies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `hamiltonian()`, `assume()`, `ppm2hz()`, `step()`, `ca_hz()`, `true()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `false()`, `eval_hz()`, `bloch_siegert()`.
