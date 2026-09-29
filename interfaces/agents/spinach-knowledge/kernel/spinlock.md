# kernel/spinlock.m

- Signature: rho = spinlock(spin_system,Lx,Ly,rho,direction)
- Source: [kernel/spinlock.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/spinlock.m)
- Wiki: [spinlock.m](https://spindynamics.org/wiki/index.php?title=spinlock.m)

## Purpose

Applies the source-described analytical approximation to spin locking: it removes spin-spin correlations and magnetisation components other than those along the selected X or Y direction. It returns the transformed state rho; it is not a time-dependent pulse simulation.

## Operation

For direction 'X', the source applies step with Ly at pi/2, homospoil with the 'destroy' option, then step with Ly at -pi/2. For direction 'Y', it uses Lx for the same three operations. direction must be 'X' or 'Y'.

## Parameters / inputs

- spin_system — Spinach system description used by step and homospoil.
- Lx, Ly — X- and Y-magnetisation operators for the spins to be locked; they must be numeric matrices of equal dimensions.
- rho — numeric state vector or bookshelf stack with the row dimension matching Lx and Ly.
- direction — 'X' or 'Y'.

## Output

- rho — the state after the selected approximation, in the representation supplied to the helper.

## Side effects

This source contains no display or file-write calls. Its stated result is the returned rho, formed through step and homospoil. No frequency or field unit is specified by this function.
