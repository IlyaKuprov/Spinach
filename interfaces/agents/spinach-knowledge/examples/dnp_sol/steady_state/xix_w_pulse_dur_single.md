# examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_w_pulse_dur_single.m`
- Signature: `xix_w_pulse_dur_single()`
- Total lines: 102

## Purpose

2D parameter scan of XiX DNP in the steady state. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: W-band magnet; implemented by `sys.magnet=3.4`.
- Lines 14-15: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 17-18: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 21-22: Spin temperature; implemented by `inter.temperature=80`.
- Lines 24-25: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 28-29: Get electron-nuclear distance; implemented by `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 38-39: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 41-42: Electron pulse duration grid, s; implemented by `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 44-45: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 53-54: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 57-58: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 60-61: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 68-69: Preallocate steady state DNP array; implemented by `dnp=zeros([numel(parameters.el_offs) numel(pulse_durs)],'like',1i)`.
- Lines 71-72: Get a figure; implemented by `kfigure()`.
- Lines 74-75: Over pulse durations; implemented by `for m=1:numel(pulse_durs)`.
- Lines 77-78: Set pulse duration; implemented by `parameters.pulse_dur=pulse_durs(m)`.

### Control flow inferred from the code

- Line 75: `for` loop over `m=1:numel(pulse_durs)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 18: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 19: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 22: computes `inter.temperature` using `inter.temperature=80`.
- Lines 25: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 29: computes `xyz` using `xyz=cell2mat(inter.coordinates); r_en=xyz(2,3)`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 39: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 42: computes `pulse_durs` using `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 45: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 46-47: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r_en,bet)`.
- Lines 48: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.
- Lines 49: computes `inter.r2_rates` using `inter.r2_rates={200e3 50e3}`.
- Lines 50: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 51: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.

## Implementation structure

- 2D parameter scan of XiX DNP in the steady state.
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
- Electron pulse duration grid, s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cell2mat()`, `xyz()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `kfigure()`, `pulse_durs()`, `dnp()`, `powder()`, `set()`, `kylabel()`, `kcolourbar()`, `kxlabel()`, `savefig()`.
