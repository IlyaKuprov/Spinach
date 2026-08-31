# examples/nmr_solids/cp_contact_mas_nhh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_contact_mas_nhh.m`
- Signature: `cp_contact_mas_nhh()`
- Total lines: 77

## Purpose

Cross-polarisation experiment in the doubly rotating frame. A single nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius sphere around it. Spinning powder simulation using a restricted Lio- uville space up to, and including, three-spin correlations. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=9.394`.
- Lines 18-20: Interactions; implemented by `inter.zeeman.scalar={0.6745 1.0368 -0.1495 -0.3171 0.7233 0.4882 3.2662 0.0317 0.0000}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Relevant operators; implemented by `Nx=operator(spin_system,'Lx','15N')`.
- Lines 53-54: MAS parameters; implemented by `parameters.rate=10000`.
- Lines 67-68: Simulation; implemented by `fid=singlerot(spin_system,@cp_contact_hard,parameters,'nmr')`.
- Lines 70-71: Plot the answer; implemented by `time_axis=[0 cumsum(parameters.time_steps)]`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 15-16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H', '1H','1H','1H','1H','15N'}`.
- Lines 19-20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.6745 1.0368 -0.1495 -0.3171 0.7233 0.4882 3.2662 0.0317 0.0000}`.
- Lines 21: computes `inter.coordinates` using `inter.coordinates={[-2.51887819 -0.99807636 0.87365165]`.
- Lines 30: computes `inter.temperature` using `inter.temperature=298`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 35: computes `bas.level` using `bas.level=3`.
- Lines 38: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 39: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 40: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 41: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `Nx` using `Nx=operator(spin_system,'Lx','15N')`.
- Lines 49: computes `Ny` using `Ny=operator(spin_system,'Ly','15N')`.
- Lines 50: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 51: computes `Hy` using `Hy=operator(spin_system,'Ly','1H')`.
- Lines 54: computes `parameters.rate` using `parameters.rate=10000`.

## Implementation structure

- Cross-polarisation experiment in the doubly rotating frame. A single
- nitrogen-15 in a bath of 8 protons scattered on a 2 Angstrom radius
- sphere around it. Spinning powder simulation using a restricted Lio-
- uville space up to, and including, three-spin correlations.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Relevant operators
- MAS parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `singlerot()`, `cumsum()`, `kfigure()`, `kylabel()`, `kxlabel()`.
