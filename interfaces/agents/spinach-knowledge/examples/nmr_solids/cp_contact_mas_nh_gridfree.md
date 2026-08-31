# examples/nmr_solids/cp_contact_mas_nh_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_contact_mas_nh_gridfree.m`
- Signature: `cp_contact_mas_nh_gridfree()`
- Total lines: 64

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 and a single proton. Spinning powder simulation starting from the thermal equilibrium using the grid-free version of the Fok- ker-Planck formalism. Calculation time: minutes with a Tesla A100 GPU, much longer otherwise.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: System specification; implemented by `sys.magnet=9.394`.
- Lines 18-19: Interactions; implemented by `inter.zeeman.scalar={0.00 0.00}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-32: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 41-42: MAS parameters; implemented by `parameters.rate=10000`.
- Lines 54-55: Simulation; implemented by `fid=gridfree(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 57-58: Plot the answer; implemented by `time_axis=[0 cumsum(parameters.time_steps)]`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'15N','1H'}`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.00 0.00}`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 22: computes `inter.temperature` using `inter.temperature=298`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 37: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 38: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 39: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 42: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 43: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 44: computes `parameters.max_rank` using `parameters.max_rank=42`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'1H','15N'}`.
- Lines 46: computes `parameters.irr_powers` using `parameters.irr_powers=[5e4*ones(1,100)`.
- Lines 48: computes `parameters.irr_opers` using `parameters.irr_opers={Hy Nx}`.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 and a single proton. Spinning powder simulation starting
- from the thermal equilibrium using the grid-free version of the Fok-
- ker-Planck formalism.
- Calculation time: minutes with a Tesla A100 GPU,
- much longer otherwise.
- System specification
- Interactions
- Basis set
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `gridfree()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
