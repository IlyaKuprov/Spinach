# etc/textbook/levelpop.m

- Signature: `[E,P,dP]=levelpop(isotope,field,temperature)`

## Purpose

Calculates Zeeman-level energies and their thermal populations for the spin specified by `isotope`, in a static magnetic field and at a specified temperature. The spin multiplicity and magnetogyric ratio are obtained from `spin`.

## Method

The function builds the Zeeman Hamiltonian `H=-mg_ratio*field*S.z`. It returns its energy levels in units of `k_B T`, then evaluates and normalizes the Boltzmann factors `exp(-E)`. The population calculation is shifted by the minimum energy before exponentiation for numerical stability. `dP` contains the signed differences `P(j)-P(j+1)` for adjacent entries in the returned population vector. Because the Zeeman Hamiltonian includes the signed magnetogyric ratio, the level ordering and population differences depend on its sign.

## Inputs

- `isotope` — character array identifying the spin (for example, `'1H'`, `'13C'`, or `'E'`).
- `field` — real scalar magnetic field in tesla.
- `temperature` — non-zero real scalar spin temperature in kelvin.

## Outputs

- `E` — vector of Zeeman energies divided by `k_B T`.
- `P` — normalized vector of level populations.
- `dP` — signed population differences for adjacent vector entries.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=levelpop.m).
