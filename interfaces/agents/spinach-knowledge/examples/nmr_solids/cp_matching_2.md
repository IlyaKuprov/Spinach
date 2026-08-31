# examples/nmr_solids/cp_matching_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_2.m`
- Signature: `cp_matching_2()`
- Total lines: 70

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus. The experiment is run with low power on the 15N nucleus, showing matching con- dition reflections with the opposite phase. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=9.394`.
- Lines 16-17: Interactions; implemented by `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 35-36: Power levels; implemented by `powers=linspace(0e3,30e3,50)`.
- Lines 38-39: Experiment parameters; implemented by `parameters.rate=10000`.
- Lines 50-51: Parallel loop over power levels; implemented by `parfor n=1:numel(powers)`.
- Lines 53-54: MAS parameters; implemented by `localpar=parameters`.
- Lines 58-59: Simulation; implemented by `fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr')`.
- Lines 64-65: Plotting; implemented by `kfigure(); plot(powers'/1e3,cp); xlim tight`.

### Control flow inferred from the code

- Line 51: `parfor` loop over `n=1:numel(powers)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','15N'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 18: computes `inter.coordinates` using `inter.coordinates={[-1.11551509 1.65289357 -1.19927242]`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 31: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 32: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 33: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 36: computes `powers` using `powers=linspace(0e3,30e3,50)`.
- Lines 39: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 40: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 41: computes `parameters.max_rank` using `parameters.max_rank=3`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'1H','15N'}`.
- Lines 43: computes `parameters.irr_opers` using `parameters.irr_opers={Hx Nx}`.
- Lines 44: computes `parameters.exc_opers` using `parameters.exc_opers={0*Hy 0*Ny}`.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus. The experiment
- is run with low power on the 15N nucleus, showing matching con-
- dition reflections with the opposite phase.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Relevant operators
- Power levels
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powers()`, `singlerot()`, `fid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
