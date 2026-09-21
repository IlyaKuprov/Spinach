# kernel/spinlock.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/spinlock.m`
- Signature: `rho=spinlock(spin_system,Lx,Ly,rho,direction)`
- Total lines: 85

## Purpose

Analytical approximation to a spin locking process. This function oblite- rates all spin-spin correlations and all magnetization components other than those along the indicated direction. Syntax: rho=spinlock(spin_system,Lx,Ly,rho,direction)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `homospoil()`, `ischar()`, `ismember()`, `ishermitian()`, `all()`.
