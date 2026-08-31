# examples/dnp_sol/steady_state/xix_w_field_profile_single.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_w_field_profile_single.m`
- Signature: `xix_w_field_profile_single()`
- Total lines: 82

## Purpose

Simulation of XiX DNP field profile in the steady state, a single spin system without ensemble averaging. Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: W-band magnet; implemented by `sys.magnet=3.4`.
- Lines 15-16: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 18-19: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 22-23: Spin temperature; implemented by `inter.temperature=80`.
- Lines 25-26: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29-30: Get electron-nuclear distance; implemented by `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 32-33: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 41-42: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 45-46: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 48-49: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 51-52: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 55-56: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 58-59: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 69-70: Run the steady state simulation; implemented by `dnp=powder(spin_system,@xixdnp_steady,parameters,'esr')`.
- Lines 72-73: Plotting; implemented by `kfigure(); plot(parameters.el_offs/1e6,real(dnp))`.
- Lines 78-79: Save the figure; implemented by `savefig(gcf,'xix_w_field_profile_single.fig')`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 19: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 20: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 23: computes `inter.temperature` using `inter.temperature=80`.
- Lines 26: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 30: computes `xyz` using `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 34-35: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r_en,bet)`.
- Lines 36: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 37: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.
- Lines 38: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 39: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 42: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 43: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 46: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 49: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 52: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Simulation of XiX DNP field profile in the steady state,
- a single spin system without ensemble averaging.
- Calculation time: seconds
- W-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Relaxation rates, distance and ori. dep. R1n
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cell2mat()`, `xyz()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `powder()`, `kfigure()`, `kylabel()`, `kxlabel()`, `savefig()`.
