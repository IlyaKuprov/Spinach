# kernel/thermalize.m

- Signature: `R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)`

## Purpose

Modifies the relaxation superoperator to drive the system to the user-specified target state using the inhomogeneous master equation (IME) formalism, or to the equilibrium state of the lab-frame Hamiltonian at the specified temperature using the DiBari-Levitt formalism.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.
- `R` - symmetric negative-definite relaxation superoperator that drives the system toward the zero state vector; it may be obtained from `relaxation.m` when `inter.equilibrium` is `'zero'`.
- `HLSPS` - lab-frame Hamiltonian left-side product superoperator, available from `hamiltonian.m` (and `orientation.m` if needed); pass an empty array for IME.
- `T` - absolute temperature; pass an empty array for IME.
- `rho_eq` - thermal equilibrium state; pass an empty array for DiBari-Levitt.
- `method` - `'IME'` for the inhomogeneous master equation or `'dibari'` for DiBari-Levitt thermalisation.

## Output

- `R` - thermalized relaxation superoperator.

## Behavior

The function validates the inputs and rejects an `R` that already acts on the unit state within its tolerance. For IME, it constructs the unit state according to the Liouville-space formalism and applies the correction `R = R - kron(U', R*rho_eq)`. For DiBari-Levitt, it computes `beta = hbar/(kbol*T)` and applies `R = R*propagator(spin_system,HLSPS,1i*beta)`.

IME requires the population of the unit state in the state vector to be exactly 1; the function cannot check or enforce this. The DiBari-Levitt method is computationally expensive, but tends to work better than IME, particularly in exotic regimes.
