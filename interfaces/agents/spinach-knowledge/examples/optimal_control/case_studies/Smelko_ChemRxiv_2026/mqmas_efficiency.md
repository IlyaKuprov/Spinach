# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m`
- Signature: `mqmas_efficiency()`
- Total lines: 141

## Purpose

Efficiency of the z-filtered 27Al MQMAS pulse sequence with hard pulses and with the optimal control pulses produced by the other examples in this folder. Reproduces, using Spinach, the sequence efficiency calculation from The nucleus has the quadrupolar coupling of aluminium acetylace- tonate (CQ=3.2 MHz, eta=0.16) and is spun at 12.5 kHz in a 400 MHz magnet; the quadrupolar interaction is taken to second order in t

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Coherence order, 3 or 5; implemented by `mq_order=5`.
- Lines 36-37: 400 MHz magnet; implemented by `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 40-41: Quadrupolar coupling and shielding anisotropy; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(3.2e6,0.16,5/2,[0 0 0])`.
- Lines 45-46: Hilbert space formalism; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 54-55: Rotor phase resolved drift Hamiltonians for the whole sequence; implemented by `parameters.spins={'27Al'}`.
- Lines 64-65: Initial state, Iz; implemented by `rho_init=state(spin_system,'Lz','27Al')`.
- Lines 68-69: Control operators; implemented by `Lx=operator(spin_system,'Lx','27Al')`.
- Lines 72-73: Hard pulse sequence, 100 kHz excitation and conversion pulses; implemented by `if mq_order==3`.
- Lines 82-83: Optimal control pulse sequence; implemented by `exc=load(['mq_exc_' num2str(mq_order) 'q.mat'],'pulse','pulse_dt')`.
- Lines 88-89: Loop over the two sequences; implemented by `sequences={hard_pulses,oc_pulses}; efficiency=zeros(1,2)`.
- Lines 92-93: Parallel loop over the ensemble; implemented by `signals=zeros(1,numel(drifts)); pulses=sequences{s}`.
- Lines 96-97: Start with the initial state at time zero; implemented by `rho=rho_init; current_time=0`.
- Lines 99-100: Loop over the three pulses; implemented by `for k=1:3`.
- Lines 102-103: Loop over the slices of the pulse; implemented by `for m=1:numel(pulses{k}{2})`.
- Lines 105-106: Drift Hamiltonian at the slice midpoint; implemented by `tick_idx=floor((current_time+pulses{k}{2}(m)/2)/tick_dt)+1`.
- Lines 108-111: Take a time step; implemented by `rho=step(spin_system,drifts{n}{tick_idx}+ pulses{k}{1}(1,m)*Lx+pulses{k}{1}(2,m)*Ly, rho,pulses{k}{2}(m))`.
- Lines 116-117: Coherence filters after excitation and conversion; implemented by `if k==1`.

### Control flow inferred from the code

- Line 73: conditional branch on `mq_order==3`.
- Line 90: `for` loop over `s=1:2`.
- Line 94: `parfor` loop over `n=1:numel(drifts)`.
- Line 100: `for` loop over `k=1:3`.
- Line 103: `for` loop over `m=1:numel(pulses{k}{2})`.
- Line 117: conditional branch on `k==1`.

### Key state/data transformations

- Lines 34: computes `mq_order` using `mq_order=5`.
- Lines 37: computes `sys.magnet` using `sys.magnet=2*pi*400e6/spin('1H')`.
- Lines 38: computes `sys.isotopes` using `sys.isotopes={'27Al'}`.
- Lines 41: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.2e6,0.16,5/2,[0 0 0])`.
- Lines 42: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-5 -5 10]}`.
- Lines 43: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 46: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 47: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 50: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 55: computes `parameters.spins` using `parameters.spins={'27Al'}`.
- Lines 56: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 57: computes `parameters.grid` using `parameters.grid='rep_2ang_400pts_sph'`.
- Lines 58: computes `parameters.n_ticks` using `parameters.n_ticks=160`.
- Lines 59: computes `parameters.n_phases` using `parameters.n_phases=32`.
- Lines 60: computes `parameters.n_slices` using `parameters.n_slices=1060`.
- Lines 61: computes `drifts` using `drifts=mqmas_drifts(spin_system,parameters)`.
- Lines 62: computes `tick_dt` using `tick_dt=0.5e-6`.
- Lines 65: computes `rho_init` using `rho_init=state(spin_system,'Lz','27Al')`.

## Implementation structure

- Efficiency of the z-filtered 27Al MQMAS pulse sequence with hard
- pulses and with the optimal control pulses produced by the other
- examples in this folder. Reproduces, using Spinach, the sequence
- efficiency calculation from
- The nucleus has the quadrupolar coupling of aluminium acetylace-
- tonate (CQ=3.2 MHz, eta=0.16) and is spun at 12.5 kHz in a 400
- MHz magnet; the quadrupolar interaction is taken to second order
- in the rotating frame. The sequence is: excitation pulse, +MQ/-MQ
- coherence filter, conversion pulse, population filter, central
- transition selective pulse, and detection of the central transi-
- tion single-quantum coherence. The efficiency is the modulus of
- the detected element of the density matrix, normalised to the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `assume()`, `mqmas_drifts()`, `state()`, `operator()`, `load()`, `num2str()`, `step()`, `coherence()`.
