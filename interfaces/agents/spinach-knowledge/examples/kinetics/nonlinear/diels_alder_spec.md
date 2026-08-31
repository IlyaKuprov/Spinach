# examples/kinetics/nonlinear/diels_alder_spec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/nonlinear/diels_alder_spec.m`
- Signature: `diels_alder_spec()`
- Total lines: 265

## Purpose

Repeated pulse-acquire experiment during the Diels-Alder cyclo- addition of acetylene to butadiene, demonstrating the non-linear kinetics module. Calculation time: hours, GPU is hard-coded.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: DFT import options; implemented by `options.min_j=2.0`.
- Lines 16-17: Load and display acetylene (substance A); implemented by `props_a=gparse('acetylene.out')`.
- Lines 25-26: Load and display butadiene (substance B); implemented by `props_b=gparse('butadiene.out')`.
- Lines 34-35: Load and display cyclohexadiene (substance C); implemented by `props_c=gparse('cyclohexadiene.out')`.
- Lines 43-44: Add natural abundance ethanol (substance D); implemented by `sys_d.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 54-56: Merge the spin systems; implemented by `[sys,inter]=merge_inp({sys_a, sys_b, sys_c, sys_d}, {inter_a,inter_b,inter_c,inter_d})`.
- Lines 58-59: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 61-62: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 64-65: Relaxation theory parameters; implemented by `inter.relaxation={'redfield','t1_t2'}`.
- Lines 77-78: Chemical parts and unit concentrations; implemented by `inter.chem.parts={1:2, 3:8, 9:16, 17:22}`.
- Lines 81-82: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 85-86: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 89-90: Reaction rate constant; implemented by `rrc=25.0`.
- Lines 92-93: 2nd order reaction generator; implemented by `K=@(t,x)(1i*[-rrc*x(2) 0 0 0`.
- Lines 98-99: Kinetic time grid, 10 seconds; implemented by `chem_nsteps=100; chem_tmax=10`.
- Lines 103-104: Preallocate concentration trajectory; implemented by `conc_traj=zeros(4,chem_nsteps+1)`.
- Lines 106-107: Initial concentrations, mol/L; implemented by `conc_traj(:,1)=[1e-2; 2e-2; 0.0; 17.1]`.
- Lines 109-110: Stage 1: concentration dynamics; implemented by `for n=1:chem_nsteps`.

### Control flow inferred from the code

- Line 110: `for` loop over `n=1:chem_nsteps`.
- Line 152: `for` loop over `n=1:chem_nsteps`.
- Line 195: `parfor` loop over `n=0:8`.
- Line 216: `for` loop over `k=1:parameters.nsteps`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=2.0`.
- Lines 14: computes `options.style` using `options.style='harmonics'`.
- Lines 17: computes `props_a` using `props_a=gparse('acetylene.out')`.
- Lines 18: computes `[sys_a,inter_a]` using `[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8,options)`.
- Lines 26: computes `props_b` using `props_b=gparse('butadiene.out')`.
- Lines 27: computes `[sys_b,inter_b]` using `[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8,options)`.
- Lines 35: computes `props_c` using `props_c=gparse('cyclohexadiene.out')`.
- Lines 36: computes `[sys_c,inter_c]` using `[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8,options)`.
- Lines 44: computes `sys_d.isotopes` using `sys_d.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 45-46: computes `inter_d.zeeman.matrix` using `inter_d.zeeman.matrix={1.26, 1.26, 1.26, 3.69, 3.69, 2.61}*eye(3)`.
- Lines 47: computes `inter_d.coordinates` using `inter_d.coordinates={[]; []; []; []; []; []}`.
- Lines 48: computes `inter_d.coupling.scalar` using `inter_d.coupling.scalar=zeros(6,6)`.
- Lines 49: computes `inter_d.coupling.scalar(1,[4 5])` using `inter_d.coupling.scalar(1,[4 5])=7.0`.
- Lines 50: computes `inter_d.coupling.scalar(2,[4 5])` using `inter_d.coupling.scalar(2,[4 5])=7.0`.
- Lines 51: computes `inter_d.coupling.scalar(3,[4 5])` using `inter_d.coupling.scalar(3,[4 5])=7.0`.
- Lines 55-56: computes `[sys,inter]` using `[sys,inter]=merge_inp({sys_a, sys_b, sys_c, sys_d}, {inter_a,inter_b,inter_c,inter_d})`.
- Lines 59: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 62: computes `sys.enable` using `sys.enable={'greedy'}`.

## Implementation structure

- Repeated pulse-acquire experiment during the Diels-Alder cyclo-
- addition of acetylene to butadiene, demonstrating the non-linear
- kinetics module.
- Calculation time: hours, GPU is hard-coded.
- DFT import options
- Load and display acetylene (substance A)
- Load and display butadiene (substance B)
- Load and display cyclohexadiene (substance C)
- Add natural abundance ethanol (substance D)
- Merge the spin systems
- Magnet field
- Greedy parallelisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `g2spinach()`, `kfigure()`, `scale_figure()`, `subplot()`, `cst_display()`, `camorbit()`, `ktitle()`, `num2cell()`, `merge_inp()`, `create()`, `basis()`, `conc_traj()`, `step()`, `kxlabel()`, `kylabel()`.
