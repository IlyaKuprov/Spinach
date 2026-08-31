# examples/dnp_sol/solid_effect_timedep_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_timedep_1.m`
- Signature: `solid_effect_timedep_1()`
- Total lines: 91

## Purpose

A simulation of solid effect DNP for a tilted linear chain of three protons positioned at distances 7, 10 and 14 Angstrom from electron located at the origin. Weizmann DNP relaxation model is used with second order Krylov-Bogolyubov average Hamiltonian theory. Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnetic field; implemented by `sys.magnet=3.4`.
- Lines 20-21: Spin system; implemented by `sys.isotopes={'E','1H','1H','1H'}`.
- Lines 28-29: Relaxation theory; implemented by `inter.relaxation={'weizmann'}`.
- Lines 46-47: Microwave power and offset; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 50-51: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 55-56: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 59-60: Experiment parameters; implemented by `parameters.theory='kb_second_order'`.
- Lines 64-65: Time domain simulation; implemented by `parameters.calc_type='time_dependence'`.
- Lines 68-70: Time domain plotting; implemented by `time_axis=linspace(0,parameters.time_step*parameters.n_steps, parameters.n_steps+1)`.
- Lines 84-85: Steady state simulation; implemented by `parameters.calc_type='steady_state'`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=2:3`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H','1H'}`.
- Lines 22: computes `R` using `R=euler2dcm(pi/6,pi/7,pi/8)`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[0 0 0 ]*R`.
- Lines 29: computes `inter.relaxation` using `inter.relaxation={'weizmann'}`.
- Lines 30: computes `inter.weiz_r1e` using `inter.weiz_r1e=1e2`.
- Lines 31: computes `inter.weiz_r1n` using `inter.weiz_r1n=0.1`.
- Lines 32: computes `inter.weiz_r2e` using `inter.weiz_r2e=1e5`.
- Lines 33: computes `inter.weiz_r2n` using `inter.weiz_r2n=1e3`.
- Lines 34: computes `inter.weiz_r1d` using `inter.weiz_r1d=zeros(4,4)`.
- Lines 35: computes `inter.weiz_r2d` using `inter.weiz_r2d=zeros(4,4)`.
- Lines 37: computes `inter.weiz_r1d(n,n+1)` using `inter.weiz_r1d(n,n+1)=0.1`.
- Lines 38: computes `inter.weiz_r1d(n+1,n)` using `inter.weiz_r1d(n+1,n)=0.1`.
- Lines 39: computes `inter.weiz_r2d(n,n+1)` using `inter.weiz_r2d(n,n+1)=0.1`.
- Lines 40: computes `inter.weiz_r2d(n+1,n)` using `inter.weiz_r2d(n+1,n)=0.1`.
- Lines 42: computes `inter.temperature` using `inter.temperature=4.2`.
- Lines 43: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 44: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.

## Implementation structure

- A simulation of solid effect DNP for a tilted linear chain of three
- protons positioned at distances 7, 10 and 14 Angstrom from electron
- located at the origin. Weizmann DNP relaxation model is used with
- second order Krylov-Bogolyubov average Hamiltonian theory.
- Calculation time: seconds
- Magnetic field
- Spin system
- Relaxation theory
- Microwave power and offset
- Basis set
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `create()`, `basis()`, `solid_effect()`, `kfigure()`, `scale_figure()`, `subplot()`, `answer()`, `kylabel()`, `kxlabel()`, `klegend()`, `set()`.
