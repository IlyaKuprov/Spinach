# examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m`
- Signature: `top_q_con_time_ensemble_r()`
- Total lines: 126

## Purpose

Simulation of TOP DNP contact time dependence in the steady state with electron-proton distance ensembles. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Q-band magnet; implemented by `sys.magnet=1.2142`.
- Lines 15-16: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 18-19: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 22-23: Spin temperature; implemented by `inter.temperature=80`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 32-33: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 35-36: Distance ensemble; implemented by `[r,wr]=gaussleg(3.5,20,3)`.
- Lines 38-39: TOP loop count; implemented by `loop_counts=1:256`.
- Lines 41-42: Preallocate equilibrium DNP value array; implemented by `dnp_a=zeros([numel(loop_counts) numel(r)],'like',1i)`.
- Lines 45-46: Over distances; implemented by `for n=1:numel(r)`.
- Lines 48-49: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 52-54: Relaxation rates, distance and orientation dependence provided using a function handle; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 62-63: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 66-67: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 69-70: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 76-77: Over loop counts; implemented by `parfor m=1:numel(loop_counts)`.
- Lines 79-80: Localise parameters; implemented by `localpar=parameters`.

### Control flow inferred from the code

- Line 46: `for` loop over `n=1:numel(r)`.
- Line 77: `parfor` loop over `m=1:numel(loop_counts)`.

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
- Lines 36: computes `[r,wr]` using `[r,wr]=gaussleg(3.5,20,3)`.
- Lines 39: computes `loop_counts` using `loop_counts=1:256`.
- Lines 42: computes `dnp_a` using `dnp_a=zeros([numel(loop_counts) numel(r)],'like',1i)`.
- Lines 43: computes `dnp_b` using `dnp_b=zeros([numel(loop_counts) numel(r)],'like',1i)`.
- Lines 49: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 54: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 55-56: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r(n),bet)`.
- Lines 57: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 58: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.

## Implementation structure

- Simulation of TOP DNP contact time dependence in the
- steady state with electron-proton distance ensembles.
- Calculation time: hours.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble
- TOP loop count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `loop_counts()`, `dnp_a()`, `powder()`, `dnp_b()`, `kfigure()`, `kylabel()`, `klegend()`, `kxlabel()`, `savefig()`.
