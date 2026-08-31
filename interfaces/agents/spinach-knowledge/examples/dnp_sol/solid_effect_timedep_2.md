# examples/dnp_sol/solid_effect_timedep_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_timedep_2.m`
- Signature: `solid_effect_timedep_2()`
- Total lines: 101

## Purpose

A simulation of solid effect DNP for a tilted linear chain of protons positioned at distances of 4+n^2 Angstrom with n=2:6 and the electron located at the origin. Weizmann DNP relaxation model is used with se- cond order Krylov-Bogolyubov average Hamiltonian theory and state space restriction to five-spin orders. Calculation time: minutes with a Tesla A100 GPU

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Magnetic field; implemented by `sys.magnet=3.4`.
- Lines 21-22: Electron; implemented by `sys.isotopes{1}='E'`.
- Lines 25-26: Nuclei; implemented by `R=euler2dcm(pi/6,pi/7,pi/8)`.
- Lines 32-33: Relaxation theory; implemented by `inter.relaxation={'weizmann'}`.
- Lines 50-51: Microwave power and offset; implemented by `parameters.mw_pwr=2*pi*250e3`.
- Lines 54-55: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 60-61: Algorithmic options; implemented by `sys.disable={'krylov'}`.
- Lines 62-65: sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 64-65: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 68-69: Experiment parameters; implemented by `parameters.theory='kb_second_order'`.
- Lines 73-74: Time domain simulation; implemented by `parameters.calc_type='time_dependence'`.
- Lines 77-79: Time domain plotting; implemented by `time_axis=linspace(0,parameters.time_step*parameters.n_steps, parameters.n_steps+1)`.
- Lines 94-95: Steady state simulation; implemented by `parameters.calc_type='steady_state'`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=2:7`.
- Line 40: `for` loop over `n=2:6`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 22: computes `sys.isotopes{1}` using `sys.isotopes{1}='E'`.
- Lines 23: computes `inter.coordinates{1}` using `inter.coordinates{1}=[0 0 0]`.
- Lines 26: computes `R` using `R=euler2dcm(pi/6,pi/7,pi/8)`.
- Lines 28: computes `sys.isotopes{n}` using `sys.isotopes{n}='1H'`.
- Lines 29: computes `inter.coordinates{n}` using `inter.coordinates{n}=[0 0 4+n^2]*R`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'weizmann'}`.
- Lines 34: computes `inter.weiz_r1e` using `inter.weiz_r1e=1e2`.
- Lines 35: computes `inter.weiz_r1n` using `inter.weiz_r1n=0.1`.
- Lines 36: computes `inter.weiz_r2e` using `inter.weiz_r2e=1e5`.
- Lines 37: computes `inter.weiz_r2n` using `inter.weiz_r2n=1e3`.
- Lines 38: computes `inter.weiz_r1d` using `inter.weiz_r1d=zeros(7,7)`.
- Lines 39: computes `inter.weiz_r2d` using `inter.weiz_r2d=zeros(7,7)`.
- Lines 41: computes `inter.weiz_r1d(n,n+1)` using `inter.weiz_r1d(n,n+1)=0.1`.
- Lines 42: computes `inter.weiz_r1d(n+1,n)` using `inter.weiz_r1d(n+1,n)=0.1`.
- Lines 43: computes `inter.weiz_r2d(n,n+1)` using `inter.weiz_r2d(n,n+1)=0.1`.
- Lines 44: computes `inter.weiz_r2d(n+1,n)` using `inter.weiz_r2d(n+1,n)=0.1`.
- Lines 46: computes `inter.temperature` using `inter.temperature=4.2`.

## Implementation structure

- A simulation of solid effect DNP for a tilted linear chain of protons
- positioned at distances of 4+n^2 Angstrom with n=2:6 and the electron
- located at the origin. Weizmann DNP relaxation model is used with se-
- cond order Krylov-Bogolyubov average Hamiltonian theory and state
- space restriction to five-spin orders.
- Calculation time: minutes with a Tesla A100 GPU
- Magnetic field
- Electron
- Nuclei
- Relaxation theory
- Microwave power and offset
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `create()`, `basis()`, `solid_effect()`, `kfigure()`, `scale_figure()`, `subplot()`, `answer()`, `kylabel()`, `kxlabel()`, `klegend()`, `set()`.
