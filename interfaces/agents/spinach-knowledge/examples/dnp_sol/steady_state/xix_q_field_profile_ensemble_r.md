# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m`
- Signature: `xix_q_field_profile_ensemble_r()`
- Total lines: 96

## Purpose

Simulation of XiX DNP field profile in the steady state with electron-proton distance ensemble averaging. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Q-band magnet; implemented by `sys.magnet=1.2142`.
- Lines 15-16: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 18-19: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 22-23: Spin temperature; implemented by `inter.temperature=80`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 32-33: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 35-36: Distance ensemble, Gauss-Legendre points; implemented by `[r,w]=gaussleg(3.5,20,3)`.
- Lines 38-39: Microwave resonance offsets, Hz; implemented by `offsets=linspace(-100e6,100e6,201)`.
- Lines 41-42: Compute DNP at each distance; implemented by `dnp=zeros([numel(offsets) numel(r)])`.
- Lines 45-46: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 49-50: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 58-59: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 62-63: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 65-66: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 75-76: Calculate shot spacing; implemented by `parameters.shot_spacing=204e-6 - 2*parameters.nloops*parameters.pulse_dur`.
- Lines 78-79: Run the steady state simulation; implemented by `dnp(:,n)=powder(spin_system,@xixdnp_steady,parameters,'esr')`.
- Lines 83-84: Integrate over the distance distribution, r^2 is the Jacobian; implemented by `dnp=sum(dnp.*reshape(r.^2,[1 numel(r)]).*reshape(w,[1 numel(w)]),2)/sum((r.^2).*w)`.

### Control flow inferred from the code

- Line 43: `for` loop over `n=1:numel(r)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=1.2142`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 19: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 20: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 23: computes `inter.temperature` using `inter.temperature=80`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 33: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 36: computes `[r,w]` using `[r,w]=gaussleg(3.5,20,3)`.
- Lines 39: computes `offsets` using `offsets=linspace(-100e6,100e6,201)`.
- Lines 42: computes `dnp` using `dnp=zeros([numel(offsets) numel(r)])`.
- Lines 46: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 50: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 51-52: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r(n),bet)`.
- Lines 53: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 54: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.
- Lines 55: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.

## Implementation structure

- Simulation of XiX DNP field profile in the steady state
- with electron-proton distance ensemble averaging.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble, Gauss-Legendre points
- Microwave resonance offsets, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `dnp()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
