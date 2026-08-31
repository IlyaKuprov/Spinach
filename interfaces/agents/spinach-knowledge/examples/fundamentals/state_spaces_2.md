# examples/fundamentals/state_spaces_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_spaces_2.m`
- Signature: `state_spaces_2()`
- Total lines: 116

## Purpose

Contributions from different orders of spin correlation to the system trajectory in the pulse-acquire 1H NMR simulation of anti-3,5-difluo- roheptane (16 spins). Different curves correspond the norms of the pro- jection of the density matrix into the subspace of one-, two-, three-, etc. spin correlations. The two traces in the lower part of the figure correspond to nine-and ten-spin correlations it is clear that for 

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet induction; implemented by `sys.magnet=11.7464`.
- Lines 19-22: Isotopes; implemented by `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 24-25: Chemical shifts; implemented by `inter.zeeman.scalar=cell(1,23)`.
- Lines 43-44: J-couplings; implemented by `inter.coupling.scalar=cell(23)`.
- Lines 80-81: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 92-93: Prevent automatic state dropout; implemented by `sys.disable={'zte'}`.
- Lines 95-99: GPU is useful here sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 98-99: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 102-103: Initial condition; implemented by `rho=state(spin_system,'L+','1H')`.
- Lines 105-106: Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 108-109: Trajectory calculation; implemented by `traj=evolution(spin_system,H,[],rho,1e-3,1000,'trajectory')`.
- Lines 111-112: Trajectory analysis by spin correlation order; implemented by `kfigure(); trajan(spin_system,traj,'correlation_order')`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=11.7464`.
- Lines 20-22: computes `sys.isotopes` using `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(1,23)`.
- Lines 26: computes `inter.zeeman.scalar{14}` using `inter.zeeman.scalar{14}= 1.0092`.
- Lines 27: computes `inter.zeeman.scalar{15}` using `inter.zeeman.scalar{15}= 1.0092`.
- Lines 28: computes `inter.zeeman.scalar{16}` using `inter.zeeman.scalar{16}= 1.0092`.
- Lines 29: computes `inter.zeeman.scalar{21}` using `inter.zeeman.scalar{21}= 1.0092`.
- Lines 30: computes `inter.zeeman.scalar{22}` using `inter.zeeman.scalar{22}= 1.0092`.
- Lines 31: computes `inter.zeeman.scalar{23}` using `inter.zeeman.scalar{23}= 1.0092`.
- Lines 32: computes `inter.zeeman.scalar{11}` using `inter.zeeman.scalar{11}= 4.6834`.
- Lines 33: computes `inter.zeeman.scalar{17}` using `inter.zeeman.scalar{17}= 4.6834`.
- Lines 34: computes `inter.zeeman.scalar{10}` using `inter.zeeman.scalar{10}= 0.0000`.
- Lines 35: computes `inter.zeeman.scalar{18}` using `inter.zeeman.scalar{18}= 0.0000`.
- Lines 36: computes `inter.zeeman.scalar{8}` using `inter.zeeman.scalar{8}= 1.7970`.
- Lines 37: computes `inter.zeeman.scalar{9}` using `inter.zeeman.scalar{9}= 1.7970`.
- Lines 38: computes `inter.zeeman.scalar{13}` using `inter.zeeman.scalar{13}= 1.6942`.
- Lines 39: computes `inter.zeeman.scalar{20}` using `inter.zeeman.scalar{20}= 1.6942`.
- Lines 40: computes `inter.zeeman.scalar{19}` using `inter.zeeman.scalar{19}= 1.6370`.

## Implementation structure

- Contributions from different orders of spin correlation to the system
- trajectory in the pulse-acquire 1H NMR simulation of anti-3,5-difluo-
- roheptane (16 spins). Different curves correspond the norms of the pro-
- jection of the density matrix into the subspace of one-, two-, three-,
- etc. spin correlations. The two traces in the lower part of the figure
- correspond to nine-and ten-spin correlations it is clear that for
- practical simulation purposes, even in the absence of relaxation, only
- correlations of up to eight spins need to be accounted for.
- Calculation time: minutes, faster with a GPU.
- Magnet induction
- Isotopes
- Chemical shifts

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `false()`, `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `evolution()`, `kfigure()`, `trajan()`, `ylim()`, `set()`.
