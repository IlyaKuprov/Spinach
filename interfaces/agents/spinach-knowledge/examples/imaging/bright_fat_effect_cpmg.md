# examples/imaging/bright_fat_effect_cpmg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/bright_fat_effect_cpmg.m`
- Signature: `bright_fat_effect_cpmg()`
- Total lines: 90

## Purpose

Bright fat effect under CPMG echo train -magnetisation losses are greater in MRI experiments on J-coupled systems because co- herences are lost in the depths of the Hilbert space. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic induction; implemented by `sys.magnet=3.0`.
- Lines 15-16: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 18-19: Chemical shifts; implemented by `inter.zeeman.scalar={1.0, 2.0, 3.0, 1.0, 2.0, 3.0}`.
- Lines 21-22: J-coupling; implemented by `inter.coupling.scalar=cell(6,6)`.
- Lines 30-32: Spins 1,2,3 are molecule A and 4,5,6 are molecule B; implemented by `inter.chem.parts={[1 2 3],[4 5 6]}`.
- Lines 34-35: Kinetic rate matrix (Hz); implemented by `inter.chem.rates=[0 0; 0 0]`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 45-49: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 52-53: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 59-60: Sample geometry; implemented by `parameters.dims=[0.30 0.25]`.
- Lines 64-65: Relaxation phantoms and operators; implemented by `parameters.rlx_ph={}; parameters.rlx_op={}`.
- Lines 67-68: Initial and detection state phantoms; implemented by `load('../../etc/phantoms/bright_fat_left.mat','left')`.
- Lines 76-77: No diffusion or flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 81-82: Run the simulation; implemented by `mri=imaging(spin_system,@cpmg_dec,parameters)`.
- Lines 84-85: Plotting; implemented by `kfigure(); surf(abs(mri)); set(gca,'XDir','reverse')`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=3.0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0, 2.0, 3.0, 1.0, 2.0, 3.0}`.
- Lines 22: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6,6)`.
- Lines 23: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=0`.
- Lines 24: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=0`.
- Lines 25: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=0`.
- Lines 26: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=11`.
- Lines 27: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=17`.
- Lines 28: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=23`.
- Lines 32: computes `inter.chem.parts` using `inter.chem.parts={[1 2 3],[4 5 6]}`.
- Lines 35: computes `inter.chem.rates` using `inter.chem.rates=[0 0; 0 0]`.
- Lines 36: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 43: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 49: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 53: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Bright fat effect under CPMG echo train -magnetisation losses
- are greater in MRI experiments on J-coupled systems because co-
- herences are lost in the depths of the Hilbert space.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin system
- Chemical shifts
- J-coupling
- Spins 1,2,3 are molecule A
- and 4,5,6 are molecule B
- Kinetic rate matrix (Hz)
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `load()`, `state()`, `imaging()`, `kfigure()`, `set()`, `kxlabel()`, `kylabel()`, `ktitle()`.
