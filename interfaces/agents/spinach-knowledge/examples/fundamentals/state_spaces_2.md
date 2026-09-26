# examples/fundamentals/state_spaces_2.m

- Signature: `state_spaces_2()`

## Purpose

Contributions from different orders of spin correlation to the system trajectory in the pulse-acquire 1H NMR simulation of anti-3,5-difluo- roheptane (16 spins). Different curves correspond the norms of the pro- jection of the density matrix into the subspace of one-, two-, three-, etc. spin correlations. The two traces in the lower part of the figure correspond to nine-and ten-spin correlations it is clear that for 

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
