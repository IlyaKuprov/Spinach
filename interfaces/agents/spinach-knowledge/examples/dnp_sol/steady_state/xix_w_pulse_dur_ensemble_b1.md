# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m`
- Signature: `xix_w_pulse_dur_ensemble_b1()`
- Total lines: 115

## Purpose

2D parameter scan of XiX DNP in the steady state with electron Rabi frequency ensemble. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: W-band magnet; implemented by `sys.magnet=3.4`.
- Lines 15-16: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 18-19: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 22-23: Spin temperature; implemented by `inter.temperature=80`.
- Lines 25-26: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29-30: Get electron-nuclear distance; implemented by `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 39-40: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 42-43: B1 ensemble; implemented by `[b1,wb1]=gaussleg(10e6,20e6,5)`.
- Lines 45-46: Electron pulse duration grid, s; implemented by `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 48-49: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 57-58: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 61-62: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 64-65: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 71-72: Preallocate steady state DNP array; implemented by `dnp=zeros([numel(parameters.el_offs) numel(pulse_durs) numel(b1)],'like',1i)`.
- Lines 74-75: Over B1 fields; implemented by `for k=1:numel(b1)`.
- Lines 77-78: Set electron nutation frequency; implemented by `parameters.irr_powers=b1(k)`.

### Control flow inferred from the code

- Line 75: `for` loop over `k=1:numel(b1)`.
- Line 81: `parfor` loop over `m=1:numel(pulse_durs)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 19: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 20: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 23: computes `inter.temperature` using `inter.temperature=80`.
- Lines 26: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 30: computes `xyz` using `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 40: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 43: computes `[b1,wb1]` using `[b1,wb1]=gaussleg(10e6,20e6,5)`.
- Lines 46: computes `pulse_durs` using `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 49: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 50-51: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r_en,bet)`.
- Lines 52: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 53: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.
- Lines 54: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.

## Implementation structure

- 2D parameter scan of XiX DNP in the steady state with
- electron Rabi frequency ensemble.
- Calculation time: hours.
- W-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Basis set
- Propagator accuracy
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cell2mat()`, `xyz()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `pulse_durs()`, `dnp()`, `powder()`, `kfigure()`, `set()`, `kylabel()`, `kxlabel()`, `kcolourbar()`, `savefig()`.
