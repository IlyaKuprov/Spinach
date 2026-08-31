# examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m`
- Signature: `top_q_con_time_ensemble_b1_r()`
- Total lines: 160

## Purpose

Simulation of TOP DNP contact time dependence in the steady state with electron-proton distance and elec- tron Rabi frequency ensembles. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Q-band magnet; implemented by `sys.magnet=1.2142`.
- Lines 16-17: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 19-20: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 23-24: Spin temperature; implemented by `inter.temperature=80`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 33-34: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 36-37: Distance ensemble; implemented by `[r,wr]=gaussleg(3.5,20,3)`.
- Lines 39-40: TOP loop count; implemented by `loop_counts=1:256`.
- Lines 42-43: B1 ensembles; implemented by `[b1a,wb1a]=gaussleg(10e6,20e6,5)`.
- Lines 46-47: Preallocate equilibrium DNP value array; implemented by `dnp_a=zeros([numel(loop_counts) numel(r) numel(b1a)],'like',1i)`.
- Lines 50-51: Over distances; implemented by `for n=1:numel(r)`.
- Lines 53-54: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 57-59: Relaxation rates, distance and orientation dependence provided using a function handle; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 67-68: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 71-72: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 74-75: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 81-82: Over B1 fields; implemented by `for k=1:numel(b1a)`.

### Control flow inferred from the code

- Line 51: `for` loop over `n=1:numel(r)`.
- Line 82: `for` loop over `k=1:numel(b1a)`.
- Line 88: `parfor` loop over `m=1:numel(loop_counts)`.
- Line 110: `for` loop over `k=1:numel(b1b)`.
- Line 116: `parfor` loop over `m=1:numel(loop_counts)`.

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
- Lines 40: computes `loop_counts` using `loop_counts=1:256`.
- Lines 43: computes `[b1a,wb1a]` using `[b1a,wb1a]=gaussleg(10e6,20e6,5)`.
- Lines 44: computes `[b1b,wb1b]` using `[b1b,wb1b]=gaussleg(25e6,35e6,5)`.
- Lines 47: computes `dnp_a` using `dnp_a=zeros([numel(loop_counts) numel(r) numel(b1a)],'like',1i)`.
- Lines 48: computes `dnp_b` using `dnp_b=zeros([numel(loop_counts) numel(r) numel(b1b)],'like',1i)`.
- Lines 54: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 59: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 60-61: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r(n),bet)`.

## Implementation structure

- Simulation of TOP DNP contact time dependence in the
- steady state with electron-proton distance and elec-
- tron Rabi frequency ensembles.
- Calculation time: hours.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `b1a()`, `loop_counts()`, `dnp_a()`, `powder()`, `b1b()`, `dnp_b()`, `kfigure()`, `kylabel()`, `klegend()`, `kxlabel()`, `savefig()`.
