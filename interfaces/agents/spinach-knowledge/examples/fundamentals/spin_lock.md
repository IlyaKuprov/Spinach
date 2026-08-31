# examples/fundamentals/spin_lock.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/spin_lock.m`
- Signature: `spin_lock()`
- Total lines: 66

## Purpose

A spin-locking experiment on a two-spin system.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Isotopes; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 10-11: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 13-14: Chemical shifts; implemented by `inter.zeeman.scalar={1.0 1.5}`.
- Lines 16-17: Scalar couplings; implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Initial state; implemented by `rho=4*state(spin_system,'Lz','all')`.
- Lines 31-32: Observable states; implemented by `coil_x1=state(spin_system,{'Lx'},{1})`.
- Lines 39-40: Pulse operators; implemented by `Lx=operator(spin_system,'Lx','all')`.
- Lines 43-44: Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 46-47: Spin-locking field of 1.5 kHz along Y; implemented by `H=H+2*pi*1.5e3*Ly`.
- Lines 49-50: Initial 90-degree pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 52-55: Time evolution; implemented by `answer=evolution(spin_system,H,[coil_x1 coil_y1 coil_z1 coil_x2 coil_y2 coil_z2], rho,1e-4,100,'multichannel')`.
- Lines 57-58: Bloch sphere plot; implemented by `[X,Y,Z]=sphere; kfigure()`.

### Key state/data transformations

- Lines 8: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 11: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0 1.5}`.
- Lines 17: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `rho` using `rho=4*state(spin_system,'Lz','all')`.
- Lines 32: computes `coil_x1` using `coil_x1=state(spin_system,{'Lx'},{1})`.
- Lines 33: computes `coil_x2` using `coil_x2=state(spin_system,{'Lx'},{2})`.
- Lines 34: computes `coil_y1` using `coil_y1=state(spin_system,{'Ly'},{1})`.
- Lines 35: computes `coil_y2` using `coil_y2=state(spin_system,{'Ly'},{2})`.
- Lines 36: computes `coil_z1` using `coil_z1=state(spin_system,{'Lz'},{1})`.
- Lines 37: computes `coil_z2` using `coil_z2=state(spin_system,{'Lz'},{2})`.
- Lines 40: computes `Lx` using `Lx=operator(spin_system,'Lx','all')`.
- Lines 41: computes `Ly` using `Ly=operator(spin_system,'Ly','all')`.
- Lines 44: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- A spin-locking experiment on a two-spin system.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Initial state
- Observable states
- Pulse operators
- Hamiltonian
- Spin-locking field of 1.5 kHz along Y

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `step()`, `evolution()`, `kfigure()`, `plot3()`, `answer()`.
