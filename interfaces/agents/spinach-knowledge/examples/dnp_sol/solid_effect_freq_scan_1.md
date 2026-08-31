# examples/dnp_sol/solid_effect_freq_scan_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_freq_scan_1.m`
- Signature: `solid_effect_freq_scan_1()`
- Total lines: 89

## Purpose

A scan through the microwave frequency range in a steady state DNP experiment for a single 15N labelled urea mole- cule at a specific orientation and a specific distance from a single electron. Laboratory frame DNP simulation is carried out with state space restriction to four-spin orders and a Weizmann DNP relaxation superoperator accounting for T1 and T2 and di- polar relaxation processes. Single crystal calculatio

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Magnetic field; implemented by `sys.magnet=3.4`.
- Lines 24-25: Spin system; implemented by `sys.isotopes={'E','15N','1H','1H','15N','1H','1H'}`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 40-41: Relaxation theory; implemented by `inter.relaxation={'weizmann'}`.
- Lines 52-53: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 56-57: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 70-71: Steady state simulation; implemented by `answer=crystal(spin_system,@dnp_freq_scan,parameters,'esr')`.
- Lines 73-74: Plotting; implemented by `kfigure(); scale_figure([1.75 1.75])`.

### Key state/data transformations

- Lines 22: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={'E','15N','1H','1H','15N','1H','1H'}`.
- Lines 26: computes `inter.coordinates` using `inter.coordinates={[ 0.00000000 0.00000000 10.14358975]`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 37: computes `bas.level` using `bas.level=4`.
- Lines 38: computes `bas.projections` using `bas.projections=[-2 -1 0 +1 +2]`.
- Lines 41: computes `inter.relaxation` using `inter.relaxation={'weizmann'}`.
- Lines 42: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 43: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 44: computes `inter.weiz_r1e` using `inter.weiz_r1e=1e2`.
- Lines 45: computes `inter.weiz_r1n` using `inter.weiz_r1n=0.1`.
- Lines 46: computes `inter.weiz_r2e` using `inter.weiz_r2e=1e5`.
- Lines 47: computes `inter.weiz_r2n` using `inter.weiz_r2n=1e3`.
- Lines 48: computes `inter.weiz_r1d` using `inter.weiz_r1d=1e-3*ones(7,7)`.
- Lines 49: computes `inter.weiz_r2d` using `inter.weiz_r2d=1e-3*ones(7,7)`.
- Lines 50: computes `inter.temperature` using `inter.temperature=4.2`.
- Lines 53: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- A scan through the microwave frequency range in a steady
- state DNP experiment for a single 15N labelled urea mole-
- cule at a specific orientation and a specific distance
- from a single electron.
- Laboratory frame DNP simulation is carried out with state
- space restriction to four-spin orders and a Weizmann DNP
- relaxation superoperator accounting for T1 and T2 and di-
- polar relaxation processes. Single crystal calculation.
- Calculation time: minutes
- Magnetic field
- Spin system
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `crystal()`, `kfigure()`, `scale_figure()`, `subplot()`, `answer()`, `kxlabel()`, `kylabel()`.
