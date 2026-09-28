# kernel/spinlock.m

- Signature: `rho=spinlock(spin_system,Lx,Ly,rho,direction)`

## Purpose

Applies an analytical approximation to spin locking: it obliterates spin-spin correlations and magnetization components other than those along the selected direction.

## Physical / mathematical content

The approximation retains only the selected X or Y magnetization. It is not a full simulation of spin locking; an explicit simulation requires RF terms in the system Hamiltonian.

## Numerical / algorithmic content

For direction `'X'`, the function applies a quarter-turn using `Ly`, destroys the other components with `homospoil(...,'destroy')`, and reverses the quarter-turn. For direction `'Y'`, it performs the corresponding steps using `Lx`. Other direction values are rejected.

## Parameters / inputs

- `spin_system` — spin system used by the propagator and homospoiler.
- `Lx` — X-magnetization operator on the spins to be locked.
- `Ly` — Y-magnetization operator on the spins to be locked.
- `rho` — state vector or bookshelf stack.
- `direction` — locking direction, `'X'` or `'Y'`.

## Outputs

- `rho` — state vector or bookshelf stack after the approximation.

## Notes

For an accurate simulation of a real spin-locking process, model it explicitly by adding RF terms to the system Hamiltonian.
