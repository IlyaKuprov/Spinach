# examples/nmr_solids/cp_matching_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_1.m`
- Signature: `cp_matching_1()`
- Total lines: 68

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus under MAS. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=9.394`.
- Lines 14-15: Interactions; implemented by `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 33-34: Power levels; implemented by `powers=linspace(20e3,80e3,120)`.
- Lines 36-37: Experiment parameters; implemented by `parameters.rate=10000`.
- Lines 48-49: Parallel loop over power levels; implemented by `parfor n=1:numel(powers)`.
- Lines 51-52: MAS parameters; implemented by `localpar=parameters`.
- Lines 56-57: Simulation; implemented by `fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr')`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot(powers'/1e3,cp); xlim tight`.

### Control flow inferred from the code

- Line 49: `parfor` loop over `n=1:numel(powers)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','15N'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[-1.11551509 1.65289357 -1.19927242]`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 29: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 30: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 31: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 34: computes `powers` using `powers=linspace(20e3,80e3,120)`.
- Lines 37: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 38: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 39: computes `parameters.max_rank` using `parameters.max_rank=3`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'1H','15N'}`.
- Lines 41: computes `parameters.irr_opers` using `parameters.irr_opers={Hx Nx}`.
- Lines 42: computes `parameters.exc_opers` using `parameters.exc_opers={0*Hy 0*Ny}`.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus under MAS.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Relevant operators
- Power levels
- Experiment parameters
- Parallel loop over power levels
- MAS parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `powers()`, `singlerot()`, `fid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
