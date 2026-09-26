# kernel/spinlock.m

- Signature: `rho=spinlock(spin_system,Lx,Ly,rho,direction)`

## Purpose

Analytical approximation to a spin locking process. This function oblite- rates all spin-spin correlations and all magnetization components other than those along the indicated direction. Syntax: rho=spinlock(spin_system,Lx,Ly,rho,direction)

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- Lx -X magnetization operator on the spins that
- should be locked
- Ly -Y magnetization operator on the spins that
- should be locked
- rho -state vector or a bookshelf stack thereof
- direction -direction in which the spins should be lo-
- cked, 'X' or 'Y'.

## Outputs

- rho -state vector or a bookshelf stack thereof
- Note: this is an approximation to what happens during a real spin locking
- process. If you need a very accurate simulation, you would need to
- model the spin locking explicitly by adding RF terms to the system
- Hamiltonian.

## Implementation structure

- Analytical approximation to a spin locking process. This function oblite-
- rates all spin-spin correlations and all magnetization components other
- than those along the indicated direction. Syntax:
- rho=spinlock(spin_system,Lx,Ly,rho,direction)
- Lx -X magnetization operator on the spins that
- should be locked
- Ly -Y magnetization operator on the spins that
- rho -state vector or a bookshelf stack thereof
- direction -direction in which the spins should be lo-
- cked, 'X' or 'Y'.
- Note: this is an approximation to what happens during a real spin locking
- process. If you need a very accurate simulation, you would need to
