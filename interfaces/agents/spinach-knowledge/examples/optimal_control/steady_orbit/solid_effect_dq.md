# examples/optimal_control/steady_orbit/solid_effect_dq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/steady_orbit/solid_effect_dq.m`
- Signature: `solid_effect_dq()`
- Total lines: 127

## Purpose

Panoramic optimisation for stroboscopic steady state DNP with the timing and power settings matching the XiX ex- eriment, but complete liberty is the choice of phase.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: W-band magnet; implemented by `sys.magnet=3.35316`.
- Lines 14-15: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 17-18: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 21-22: Spin temperature; implemented by `inter.temperature=80`.
- Lines 24-25: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 28-30: Get electron-nuclear distance; implemented by `r_en=norm(inter.coordinates{2}- inter.coordinates{1},2)`.
- Lines 32-33: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 41-42: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 45-46: Parallelisation settings; implemented by `sys.parallel={'processes',240}`.
- Lines 48-49: Very tight numerics; implemented by `sys.tols.prop_chop=1e-14`.
- Lines 52-53: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 56-57: Get electron control operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 60-61: Get electron offset operator; implemented by `Ez=operator(spin_system,'Lz','E')`.
- Lines 63-65: Initial state will be ignored by StStSt, but set to thermodynamic equilibrium for when the option is off; implemented by `H=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 68-70: Target state is nuclear magnetisation, relative to the thermal equilibrium; implemented by `rho_targ=state(spin_system,'Lz','1H')`.
- Lines 74-75: Powder context parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 78-79: Transmitter is at 94.0 GHz precisely; implemented by `parameters.offset=[-spin('E')*sys.magnet/(2*pi)-94.0e9, 0]`.
- Lines 81-82: Get drift Liouvillians at all orientations; implemented by `control.drifts=drifts(spin_system,@powder,parameters,'esr')`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=3.35316`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 18: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 19: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 22: computes `inter.temperature` using `inter.temperature=80`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29-30: computes `r_en` using `r_en=norm(inter.coordinates{2}- inter.coordinates{1},2)`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 34-35: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1.0e-3,52.0,r_en,bet)`.
- Lines 36: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 37: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.
- Lines 38: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 39: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 42: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 43: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 46: computes `sys.parallel` using `sys.parallel={'processes',240}`.
- Lines 49: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-14`.
- Lines 50: computes `sys.tols.stst_tol` using `sys.tols.stst_tol=1e-10`.

## Implementation structure

- Panoramic optimisation for stroboscopic steady state DNP
- with the timing and power settings matching the XiX ex-
- eriment, but complete liberty is the choice of phase.
- W-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Relaxation rates, distance and ori. dep. R1n
- Basis set
- Parallelisation settings

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `r1n_dnp()`, `create()`, `basis()`, `operator()`, `hamiltonian()`, `assume()`, `equilibrium()`, `state()`, `spin()`, `drifts()`, `false()`, `true()`, `load()`, `firf()`, `optimcon()`, `wrapTo2Pi()`.
