# examples/nmr_solids/cp_matching_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_matching_3.m`
- Signature: `cp_matching_3()`
- Total lines: 81

## Purpose

Hartmann-Hahn matching condition test for a cross-polarisation experiment between a proton and a 15N nucleus. A 2D scan of power levels at a specific spinning rate. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=9.394`.
- Lines 15-16: Interactions; implemented by `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 39-40: Power levels; implemented by `powers=linspace(0e3,50e3,50)`.
- Lines 42-43: Experiment parameters; implemented by `parameters.rate=10000`.
- Lines 54-55: Parallel loop over power levels; implemented by `cp=zeros(numel(powers),numel(powers)); kfigure()`.
- Lines 60-61: MAS parameters; implemented by `localpar=parameters`.
- Lines 65-66: Simulation; implemented by `fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr')`.
- Lines 71-73: Plot as an image; implemented by `imagesc([min(powers) max(powers)], [min(powers) max(powers)],cp)`.

### Control flow inferred from the code

- Line 56: `for` loop over `n=1:numel(powers)`.
- Line 58: `parfor` loop over `k=1:numel(powers)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','15N'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.1495 0.0000}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[-1.11551509 1.65289357 -1.19927242]`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 26: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 27: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 35: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 36: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 37: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 40: computes `powers` using `powers=linspace(0e3,50e3,50)`.
- Lines 43: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 44: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 45: computes `parameters.max_rank` using `parameters.max_rank=3`.

## Implementation structure

- Hartmann-Hahn matching condition test for a cross-polarisation
- experiment between a proton and a 15N nucleus. A 2D scan of
- power levels at a specific spinning rate.
- Calculation time: hours
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Relevant operators
- Power levels
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `kfigure()`, `powers()`, `singlerot()`, `fid()`, `kxlabel()`, `kylabel()`, `set()`.
