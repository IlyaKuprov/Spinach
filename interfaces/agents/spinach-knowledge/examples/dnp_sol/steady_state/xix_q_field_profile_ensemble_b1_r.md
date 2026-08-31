# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m`
- Signature: `xix_q_field_profile_ensemble_b1_r()`
- Total lines: 111

## Purpose

Simulation of XiX DNP field profile in the steady state with averaging over electron-proton distance and electron Rabi frequency ensemble. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Q-band magnet; implemented by `sys.magnet=1.2142`.
- Lines 16-17: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 19-20: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 23-24: Spin temperature; implemented by `inter.temperature=80`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 33-34: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 36-37: Distance and B1 ensemble, Gauss-Legendre points; implemented by `[r,wr]=gaussleg(3.5,20,3)`.
- Lines 40-41: Microwave resonance offsets, Hz; implemented by `offsets=linspace(-100e6,100e6,201)`.
- Lines 43-44: Preallocate equilibrium DNP value array; implemented by `dnp=zeros([numel(offsets) numel(r) numel(b1)],'like',1i)`.
- Lines 46-47: Over distances; implemented by `for n=1:numel(r)`.
- Lines 49-50: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 53-55: Relaxation rates, distance and orientation dependence provided using a function handle; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 63-64: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 67-68: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 70-71: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 79-80: Calculate shot spacing; implemented by `parameters.shot_spacing=204e-6 - 2*parameters.nloops*parameters.pulse_dur`.
- Lines 82-83: Over B1 fields; implemented by `for k=1:numel(b1)`.

### Control flow inferred from the code

- Line 47: `for` loop over `n=1:numel(r)`.
- Line 83: `for` loop over `k=1:numel(b1)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=1.2142`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 20: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 21: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 24: computes `inter.temperature` using `inter.temperature=80`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 34: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 37: computes `[r,wr]` using `[r,wr]=gaussleg(3.5,20,3)`.
- Lines 38: computes `[b1,wb1]` using `[b1,wb1]=gaussleg(10e6,20e6,5)`.
- Lines 41: computes `offsets` using `offsets=linspace(-100e6,100e6,201)`.
- Lines 44: computes `dnp` using `dnp=zeros([numel(offsets) numel(r) numel(b1)],'like',1i)`.
- Lines 50: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 55: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 56-57: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r(n),bet)`.
- Lines 58: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 59: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.

## Implementation structure

- Simulation of XiX DNP field profile in the steady state with
- averaging over electron-proton distance and electron Rabi
- frequency ensemble.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance and B1 ensemble, Gauss-Legendre points

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `dnp()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
