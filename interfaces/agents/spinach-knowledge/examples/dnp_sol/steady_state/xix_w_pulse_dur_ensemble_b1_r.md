# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m`
- Signature: `xix_w_pulse_dur_ensemble_b1_r()`
- Total lines: 126

## Purpose

2D parameter scan of XiX DNP in the steady state with electron-proton distance and electron Rabi frequency ensembles. Calculation time: hours.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: W-band magnet; implemented by `sys.magnet=3.4`.
- Lines 16-17: Electron and proton; implemented by `sys.isotopes={'E','1H'}`.
- Lines 19-20: Zeeman interactions (g-tensor for trityl, ppm guess for 1H); implemented by `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 23-24: Spin temperature; implemented by `inter.temperature=80`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Propagator accuracy; implemented by `sys.tols.prop_chop=1e-12`.
- Lines 33-34: Algorithmic options; implemented by `sys.disable={'hygiene'}`.
- Lines 36-37: Distance and B1 ensemble; implemented by `[r,w]=gaussleg(3.5,20,3)`.
- Lines 40-41: Electron pulse duration grid, s; implemented by `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 43-44: Microwave resonance offsets, Hz; implemented by `offsets=linspace(-230e6,205e6,101)`.
- Lines 46-48: Preallocate steady state DNP array; implemented by `dnp=zeros([numel(offsets) numel(pulse_durs) numel(r) numel(b1)],'like',1i)`.
- Lines 50-51: Over distances; implemented by `for n=1:numel(r)`.
- Lines 53-54: Cartesian coordinates; implemented by `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 57-58: Relaxation rates, distance and ori. dep. R1n; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 66-67: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 70-71: Detect the proton; implemented by `parameters.coil=state(spin_system,'Lz','1H')`.
- Lines 73-74: Experiment parameters; implemented by `parameters.spins={'E','1H'}`.
- Lines 80-81: Over B1 fields; implemented by `for k=1:numel(b1)`.

### Control flow inferred from the code

- Line 51: `for` loop over `n=1:numel(r)`.
- Line 81: `for` loop over `k=1:numel(b1)`.
- Line 87: `parfor` loop over `m=1:numel(pulse_durs)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 20: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5]}`.
- Lines 21: computes `inter.zeeman.euler` using `inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10]}`.
- Lines 24: computes `inter.temperature` using `inter.temperature=80`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `sys.tols.prop_chop` using `sys.tols.prop_chop=1e-12`.
- Lines 34: computes `sys.disable` using `sys.disable={'hygiene'}`.
- Lines 37: computes `[r,w]` using `[r,w]=gaussleg(3.5,20,3)`.
- Lines 38: computes `[b1,wb1]` using `[b1,wb1]=gaussleg(10e6,20e6,5)`.
- Lines 41: computes `pulse_durs` using `pulse_durs=linspace(2e-9,21e-9,200)`.
- Lines 44: computes `offsets` using `offsets=linspace(-230e6,205e6,101)`.
- Lines 47-48: computes `dnp` using `dnp=zeros([numel(offsets) numel(pulse_durs) numel(r) numel(b1)],'like',1i)`.
- Lines 54: computes `inter.coordinates` using `inter.coordinates={[0.000 0.000 0.000]`.
- Lines 58: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 59-60: computes `r1n_rate` using `r1n_rate=@(alp,bet,gam)r1n_dnp(sys.magnet,inter.temperature, 2.00230,1e-3,52,r(n),bet)`.
- Lines 61: computes `inter.r1_rates` using `inter.r1_rates={1e3 r1n_rate}`.

## Implementation structure

- 2D parameter scan of XiX DNP in the steady state with
- electron-proton distance and electron Rabi frequency
- ensembles.
- Calculation time: hours.
- W-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance and B1 ensemble

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `pulse_durs()`, `dnp()`, `powder()`, `kfigure()`, `set()`, `kylabel()`, `kxlabel()`, `kcolourbar()`, `savefig()`.
