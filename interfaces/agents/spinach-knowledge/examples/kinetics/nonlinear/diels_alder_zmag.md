# examples/kinetics/nonlinear/diels_alder_zmag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/nonlinear/diels_alder_zmag.m`
- Signature: `diels_alder_zmag()`
- Total lines: 170

## Purpose

Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition of acetylene to butadiene, demonstrating the non-linear kinetics module. Calculation time: minutes.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: DFT import options; implemented by `options.min_j=2.0`.
- Lines 15-16: Load and display acetylene (substance A); implemented by `props_a=gparse('acetylene.out')`.
- Lines 24-25: Load and display butadiene (substance B); implemented by `props_b=gparse('butadiene.out')`.
- Lines 33-34: Load and display cyclohexadiene (substance C); implemented by `props_c=gparse('cyclohexadiene.out')`.
- Lines 42-43: Add natural abundance ethanol (substance D); implemented by `sys_d.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 53-55: Merge the spin systems; implemented by `[sys,inter]=merge_inp({sys_a, sys_b, sys_c, sys_d}, {inter_a,inter_b,inter_c,inter_d})`.
- Lines 57-58: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 60-61: Chemical parts and unit concentrations; implemented by `inter.chem.parts={1:2, 3:8, 9:16, 17:22}`.
- Lines 64-65: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 68-69: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 72-73: Initial concentrations, mol/L; implemented by `A0=1e-2; B0=2e-2; C0=0; D0=0.1`.
- Lines 75-76: Reaction rate constant; implemented by `rrc=25.0`.
- Lines 78-79: 2nd order reaction generator; implemented by `K=@(t,x)(1i*[-rrc*x(2) 0 0 0`.
- Lines 84-85: Time grid (ten seconds); implemented by `nsteps=100; tmax=10.0; dt=tmax/nsteps`.
- Lines 88-89: Preallocate trajectory; implemented by `x=zeros(4,nsteps+1)`.
- Lines 91-92: Define initial concentrations; implemented by `x(:,1)=[A0 B0 C0 D0]'`.
- Lines 94-95: Run Lie group solver; implemented by `for n=1:nsteps`.
- Lines 99-100: Interpolate concentrations as functions of time; implemented by `A=griddedInterpolant(time_axis,x(1,:),'makima','none')`.

### Control flow inferred from the code

- Line 95: `for` loop over `n=1:nsteps`.
- Line 135: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 12: computes `options.min_j` using `options.min_j=2.0`.
- Lines 13: computes `options.style` using `options.style='harmonics'`.
- Lines 16: computes `props_a` using `props_a=gparse('acetylene.out')`.
- Lines 17: computes `[sys_a,inter_a]` using `[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8,options)`.
- Lines 25: computes `props_b` using `props_b=gparse('butadiene.out')`.
- Lines 26: computes `[sys_b,inter_b]` using `[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8,options)`.
- Lines 34: computes `props_c` using `props_c=gparse('cyclohexadiene.out')`.
- Lines 35: computes `[sys_c,inter_c]` using `[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8,options)`.
- Lines 43: computes `sys_d.isotopes` using `sys_d.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 44-45: computes `inter_d.zeeman.matrix` using `inter_d.zeeman.matrix={1.26, 1.26, 1.26, 3.69, 3.69, 2.61}*eye(3)`.
- Lines 46: computes `inter_d.coordinates` using `inter_d.coordinates={[]; []; []; []; []; []}`.
- Lines 47: computes `inter_d.coupling.scalar` using `inter_d.coupling.scalar=zeros(6,6)`.
- Lines 48: computes `inter_d.coupling.scalar(1,[4 5])` using `inter_d.coupling.scalar(1,[4 5])=7.0`.
- Lines 49: computes `inter_d.coupling.scalar(2,[4 5])` using `inter_d.coupling.scalar(2,[4 5])=7.0`.
- Lines 50: computes `inter_d.coupling.scalar(3,[4 5])` using `inter_d.coupling.scalar(3,[4 5])=7.0`.
- Lines 54-55: computes `[sys,inter]` using `[sys,inter]=merge_inp({sys_a, sys_b, sys_c, sys_d}, {inter_a,inter_b,inter_c,inter_d})`.
- Lines 58: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 61: computes `inter.chem.parts` using `inter.chem.parts={1:2, 3:8, 9:16, 17:22}`.

## Implementation structure

- Time-domain Z magnetisation dynamics in the Diels-Alder cycloaddition
- of acetylene to butadiene, demonstrating the non-linear kinetics module.
- Calculation time: minutes.
- DFT import options
- Load and display acetylene (substance A)
- Load and display butadiene (substance B)
- Load and display cyclohexadiene (substance C)
- Add natural abundance ethanol (substance D)
- Merge the spin systems
- Magnet field
- Chemical parts and unit concentrations
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `g2spinach()`, `kfigure()`, `scale_figure()`, `subplot()`, `cst_display()`, `camorbit()`, `ktitle()`, `num2cell()`, `merge_inp()`, `create()`, `basis()`, `step()`, `griddedInterpolant()`, `kxlabel()`, `kylabel()`.
