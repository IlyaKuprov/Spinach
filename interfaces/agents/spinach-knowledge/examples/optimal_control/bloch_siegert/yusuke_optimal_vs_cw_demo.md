# examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m`
- Signature: `yusuke_optimal_vs_cw_demo()`
- Total lines: 240

## Purpose

Bloch-Siegert-aware phase optimisation compared to a simple constant- phase low-power cycle. This is the control-side companion to yusuke_14n_broadening_demo.m: the task is a reduced-model surrogate for low-power offset-tolerant decoupling, formulated as an identity cycle that should preserve magnetisation across offset and B1 distributions. The "non-optimal pulse" is a constant-phase X pulse with the same RF amplitu

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `offset_profile()`, `b1_profile()`, `ensemble_score()`, `single_score()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Magnetic field corresponding to 800 MHz 1H; implemented by `sys.magnet=18.8`.
- Lines 22-23: Single-spin surrogate model; implemented by `sys.isotopes={'13C'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: States representing the identity map; implemented by `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 39-40: Control operators and drift; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 45-52: Practical low-power regime inspired by the 14N decoupling papers: 70 kHz MAS, tp ~ 10 us, and nu_14N in the 15-23 kHz range. In this reduced model, we use a 20 kHz carrier and a 10-step cycle of 10 us elements, giving a 4pi constant-phase reference pulse. The effective offset window below is intentionally narrower than the full laboratory- frame 14N dispersion because this script is a reduced control surrogate, not a full QJF/MAS quadrupolar simulation.; implemented by `rf_nominal_hz=20e3`.
- Lines 60-61: Define the phase-only optimisation problem; implemented by `control.isotopes={'13C'}`.
- Lines 76-77: Enable Bloch-Siegert corrections in the optimiser and simulator; implemented by `control.bsiegert=true()`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system_bs=optimcon(spin_system,control)`.
- Lines 82-83: Optimise a phase-modulated identity cycle; implemented by `guess=2*pi*rand(1,nsteps)`.
- Lines 86-87: Constant-phase baseline with the same duration and RF power; implemented by `phi_cw=zeros(1,nsteps)`.
- Lines 89-92: Evaluate both waveforms on a dense offset and B1 grid; implemented by `[profile_opt,mean_offset_opt]=offset_profile(spin_system_bs,H,{Lx,Ly}, control.pulse_dt,Sx,Sy,Sz,Lz,phi_opt,rf_nominal_hz,eval_offsets_hz, eval_b1_scales)`.
- Lines 106-108: Training-set figures of merit; implemented by `train_opt=ensemble_score(spin_system_bs,H,{Lx,Ly},control.pulse_dt, Sx,Sy,Sz,Lz,phi_opt,rf_nominal_hz,train_offsets_hz,train_b1_scales)`.
- Lines 113-114: Plot the result; implemented by `figure`.
- Lines 144-145: Secondary robustness figures; implemented by `figure`.
- Lines 163-164: Console summary; implemented by `fprintf('\n')`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=18.8`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'13C'}`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `Sx` using `Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2)`.
- Lines 36: computes `Sy` using `Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2)`.
- Lines 37: computes `Sz` using `Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2)`.
- Lines 40: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 41: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 42: computes `Lz` using `Lz=operator(spin_system,'Lz','13C')`.
- Lines 43: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 52: computes `rf_nominal_hz` using `rf_nominal_hz=20e3`.
- Lines 53: computes `nsteps` using `nsteps=10`.
- Lines 54: computes `pulse_dt` using `pulse_dt=10e-6*ones(1,nsteps)`.
- Lines 55: computes `train_offsets_hz` using `train_offsets_hz=linspace(-12e3,12e3,7)`.
- Lines 56: computes `train_b1_scales` using `train_b1_scales=[0.95 1.00 1.05]`.

### Local helper functions

- Line 177: `offset_profile()` — `function [score_map,mean_profile]=offset_profile(spin_system,H,controls,`. Fidelity map across offsets and B1 scales
  - Representative operation: `pulse_dt,Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)`.
  - Representative operation: `score_map=zeros(numel(b1_scales),numel(offsets_hz))`.
- Line 192: `b1_profile()` — `function [score_map,mean_profile]=b1_profile(spin_system,H,controls,`. Fidelity map across offsets and B1 scales
  - Representative operation: `pulse_dt,Sx,Sy,Sz,Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)`.
  - Representative operation: `score_map=zeros(numel(b1_scales),numel(offsets_hz))`.
- Line 207: `ensemble_score()` — `function score=ensemble_score(spin_system,H,controls,pulse_dt,Sx,Sy,Sz,`. Training-ensemble fidelity matrix
  - Representative operation: `Lz,phi_profile,rf_nominal_hz,offsets_hz,b1_scales)`.
  - Representative operation: `score=zeros(numel(b1_scales),numel(offsets_hz))`.
- Line 221: `single_score()` — `function score=single_score(spin_system,H,controls,pulse_dt,Sx,Sy,Sz,`. Build the physical waveform
  - Representative operation: `Lz,phi_profile,rf_nominal_hz,offset_hz,b1_scale)`.
  - Representative operation: `amp_profile=2*pi*rf_nominal_hz*b1_scale*ones(size(phi_profile))`.

## Implementation structure

- Bloch-Siegert-aware phase optimisation compared to a simple constant-
- phase low-power cycle. This is the control-side companion to
- yusuke_14n_broadening_demo.m: the task is a reduced-model surrogate for
- low-power offset-tolerant decoupling, formulated as an identity cycle
- that should preserve magnetisation across offset and B1 distributions.
- The "non-optimal pulse" is a constant-phase X pulse with the same RF
- amplitude and total duration. It behaves like a simple CW-style cycle:
- acceptable near the design point, but poor once offset and B1 errors
- are included. The "optimal pulse" is a phase-modulated waveform
- optimised with Bloch-Siegert corrections enabled.
- Magnetic field corresponding to 800 MHz 1H
- Single-spin surrogate model

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rng()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `true()`, `optimcon()`, `fmaxnewton()`, `offset_profile()`, `b1_profile()`, `ensemble_score()`, `cumsum()`, `subplot()`, `unwrap()`.
