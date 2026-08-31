# examples/fundamentals/derivative_tests/difdiff_bs_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/difdiff_bs_rect.m`
- Signature: `difdiff_bs_rect()`
- Total lines: 151

## Purpose

Directional derivative test for Cartesian GRAPE with Bloch-Siegert corrections.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Set the magnetic field; implemented by `sys.magnet=10.2`.
- Lines 11-12: Set isotopes; implemented by `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 14-15: Set interactions; implemented by `inter.zeeman.scalar={1.5,2.0,30.0,40.0}`.
- Lines 22-23: Set basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Build and normalise initial state; implemented by `rho_init=singlet(spin_system,1,2)`.
- Lines 34-35: Build and normalise target state; implemented by `rho_targ=state(spin_system,{'Lz'},{4})`.
- Lines 38-39: Get control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 44-45: Get offset and shift operators; implemented by `LzH=operator(spin_system,'Lz','1H')`.
- Lines 48-49: Build drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 51-52: Define control structure; implemented by `control.drifts={{H}}`.
- Lines 62-63: Define transmitter offsets for both isotopes; implemented by `control.offsets={1050,5285}`.
- Lines 66-67: Enable Bloch-Siegert corrections; implemented by `control.bsiegert=true()`.
- Lines 71-72: Run Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 74-75: Build random guess and finite difference increment; implemented by `guess=randn(numel(control.operators),numel(control.pulse_dt))/10`.
- Lines 79-80: Request analytical gradient; implemented by `[~,~,grad_anl]=grape_xy(guess,spin_system)`.
- Lines 83-84: Left waveform edge; implemented by `wave_forw=guess`.
- Lines 101-102: Right waveform edge; implemented by `wave_forw=guess`.

### Control flow inferred from the code

- Line 92: conditional branch on `rel_err<5e-6`.
- Line 110: conditional branch on `rel_err<5e-6`.
- Line 133: conditional branch on `rel_err<5e-6`.
- Line 143: conditional branch on `fail_count==0`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=10.2`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5,2.0,30.0,40.0}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 18: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=150`.
- Lines 19: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=150`.
- Lines 20: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=50`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `rho_init` using `rho_init=singlet(spin_system,1,2)`.
- Lines 35: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{4})`.
- Lines 39: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 40: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 41: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 42: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 45: computes `LzH` using `LzH=operator(spin_system,'Lz','1H')`.

## Implementation structure

- Directional derivative test for Cartesian GRAPE with Bloch-Siegert
- corrections.
- Set the magnetic field
- Set isotopes
- Set interactions
- Set basis
- Run Spinach housekeeping
- Build and normalise initial state
- Build and normalise target state
- Get control operators
- Get offset and shift operators
- Build drift Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `singlet()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `true()`, `optimcon()`, `grape_xy()`, `squeeze()`, `grad_anl()`, `wave_forw()`, `wave_back()`, `fid_forw()`, `fid_back()`.
