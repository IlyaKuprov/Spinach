# interfaces/gaussian/brokensymm.m

- Signature: `J=brokensymm(props_sing,props_trip)`

## Purpose

Estimates the exchange coupling from singlet and triplet DFT results using the Yamaguchi equation. The Hamiltonian convention is `H=-2J*(Sa.Sb)`; the routine returns a rough order-of-magnitude estimate, not a high-precision coupling.

## Method

Using the singlet and triplet energies and their squared-spin expectation values, the function evaluates `J=(E_trip-E_sing)/(S2_sing-S2_trip)`, then converts Hartree to Hz with the factor `6.57968974479e15`. This is Eq. 6 of the cited paper.

## Parameters / inputs

- `props_sing`: `gparse` output for the singlet state of the biradical; must contain `energy` and `s_sq`.
- `props_trip`: `gparse` output for the triplet state of the biradical; must contain `energy` and `s_sq`.

## Output

- `J`: rough estimate of the exchange coupling in Hz under the stated Hamiltonian convention.

## References

- Yamaguchi equation, Eq. 6: [doi:10.1063/1.5144696](https://doi.org/10.1063/1.5144696).
- [Spinach Wiki: brokensymm.m](https://spindynamics.org/wiki/index.php?title=brokensymm.m)
