# examples/dnp_sol/cross_effect_freq_scan_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/cross_effect_freq_scan_1.m`
- Signature: `cross_effect_freq_scan_1()`
- Total lines: 72

## Purpose

A simple TOTAPOL based Cross Effect DNP system. Set to repro- duce Figure 2c from Intensity differences are due to a different relaxation model and minor inconsistencies between the stated geometry and the interaction amplitudes used in the original paper. Electron rotating frame simulation using Nottingham DNP rela- xation theory detailed in Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Magnetic field; implemented by `sys.magnet=3.4`.
- Lines 26-27: Spin system; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Relaxation theory; implemented by `inter.relaxation={'nottingham'}`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 63-64: Steady state simulation; implemented by `answer=crystal(spin_system,@dnp_freq_scan,parameters,'esr')`.
- Lines 66-67: Plotting; implemented by `kfigure(); plot(linspace(-350,350,5e4),real(answer)); kgrid`.

### Key state/data transformations

- Lines 24: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023193 2.0021091 0.0000000}`.
- Lines 29: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `inter.relaxation` using `inter.relaxation={'nottingham'}`.
- Lines 39: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 40: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 41: computes `inter.nott_r1e` using `inter.nott_r1e=1e2`.
- Lines 42: computes `inter.nott_r2e` using `inter.nott_r2e=1e5`.
- Lines 43: computes `inter.nott_r1n` using `inter.nott_r1n=0.1`.
- Lines 44: computes `inter.nott_r2n` using `inter.nott_r2n=1e3`.
- Lines 45: computes `inter.temperature` using `inter.temperature=10`.
- Lines 48: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 52: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 53: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*100e3`.
- Lines 54: computes `parameters.mw_frq` using `parameters.mw_frq=2*pi*linspace(-350,350,5e4)*1e6`.

## Implementation structure

- A simple TOTAPOL based Cross Effect DNP system. Set to repro-
- duce Figure 2c from
- Intensity differences are due to a different relaxation model
- and minor inconsistencies between the stated geometry and the
- interaction amplitudes used in the original paper.
- Electron rotating frame simulation using Nottingham DNP rela-
- xation theory detailed in
- Calculation time: seconds
- Magnetic field
- Spin system
- Basis set
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `crystal()`, `kfigure()`, `kxlabel()`, `kylabel()`.
