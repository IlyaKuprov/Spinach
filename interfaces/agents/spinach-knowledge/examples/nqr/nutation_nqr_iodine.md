# examples/nqr/nutation_nqr_iodine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nqr/nutation_nqr_iodine.m`
- Signature: `nutation_nqr_iodine()`
- Total lines: 69

## Purpose

Powder NQR nutation curve for a system with a single 127I nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=0`.
- Lines 15-16: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 19-20: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Experiment parameters; implemented by `parameters.spins={'127I'}`.
- Lines 42-43: Get a figure started; implemented by `kfigure(); scale_figure([2 1]); drawnow`.
- Lines 45-46: Loop over the pulse durations; implemented by `for n=1:10`.
- Lines 48-49: Set pulse duration; implemented by `params=parameters; params.rf_dur=5e-7*n`.
- Lines 51-52: Run the simulation; implemented by `spectrum=powder(spin_system,@nqr_pa,params,'labframe')`.
- Lines 54-55: Demodulate the transmitter offset; implemented by `spectrum=exp(-1i*(2*pi*parameters.rf_frq)*params.rf_dur)*spectrum`.
- Lines 57-60: Plotting; implemented by `frq_axis=linspace(parameters.sweep(1), parameters.sweep(2), parameters.npoints)`.

### Control flow inferred from the code

- Line 46: `for` loop over `n=1:10`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'127I'}`.
- Lines 13: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(560e6,0.01,5/2,[0 0 0])`.
- Lines 16: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 17: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 21: computes `inter.damp_rate` using `inter.damp_rate=1e5`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.temperature` using `inter.temperature=298`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'127I'}`.
- Lines 32: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=[83.5e6 84.5e6]`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 35: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 36: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','127I')`.
- Lines 37: computes `parameters.Lx` using `parameters.Lx=operator(spin_system,'Lx','127I')`.

## Implementation structure

- Powder NQR nutation curve for a system with a
- single 127I nucleus.
- Calculation time: seconds
- System specification
- Formalism and basis
- Relaxation theory
- Spinach housekeeping
- Experiment parameters
- Get a figure started
- Loop over the pulse durations
- Set pulse duration
- Run the simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `operator()`, `kfigure()`, `scale_figure()`, `powder()`, `subplot()`, `ktitle()`, `num2str()`, `kxlabel()`, `ylim()`, `set()`.
