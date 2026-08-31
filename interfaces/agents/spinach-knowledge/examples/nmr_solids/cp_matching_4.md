# examples/nmr_solids/cp_matching_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_4.m`
- Signature: `cp_matching_4()`
- Total lines: 77

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus in the presence of conformational exchange between two geometries that differ by 90 degrees in the N-H vector direction. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=9.394`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.scalar={0.0 0.0 0.0 0.0}`.
- Lines 22-23: Chemical exchange; implemented by `inter.chem.parts={[1 2],[3 4]}`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 42-43: Power levels; implemented by `powers=linspace(20e3,80e3,60)`.
- Lines 45-46: Experiment parameters; implemented by `parameters.rate=10000`.
- Lines 57-58: Parallel loop over power levels; implemented by `parfor n=1:60`.
- Lines 60-61: MAS parameters; implemented by `localpar=parameters`.
- Lines 65-66: Simulation; implemented by `fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr')`.
- Lines 71-72: Plotting; implemented by `kfigure(); plot(powers'/1e3,cp)`.

### Control flow inferred from the code

- Line 58: `parfor` loop over `n=1:60`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','15N','1H','15N'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 0.0 0.0 0.0}`.
- Lines 19: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 0 2]`.
- Lines 23: computes `inter.chem.parts` using `inter.chem.parts={[1 2],[3 4]}`.
- Lines 24: computes `inter.chem.rates` using `inter.chem.rates=[-5000 +5000`.
- Lines 26: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 38: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 39: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 40: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 43: computes `powers` using `powers=linspace(20e3,80e3,60)`.
- Lines 46: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 47: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 48: computes `parameters.max_rank` using `parameters.max_rank=3`.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus in the presence of
- conformational exchange between two geometries that differ by
- 90 degrees in the N-H vector direction.
- Calculation time: seconds
- System specification
- Interactions
- Chemical exchange
- Basis set
- Spinach housekeeping
- Relevant operators
- Power levels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powers()`, `singlerot()`, `fid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
