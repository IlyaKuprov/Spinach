# etc/textbook/levelpop.m

- MATLAB implementation: [etc/textbook/levelpop.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/levelpop.m)

## Purpose

Computes the energy levels and equilibrium populations for one isotope in a static field at a specified spin temperature. This is a single-spin Zeeman calculation; it does not assemble couplings or a multi-spin Hamiltonian.

## Use

```matlab
[E,P,dP]=levelpop(isotope,field,temperature)
```

All arguments are required; no defaults are defined.

- `isotope` — character array naming a Spinach isotope; source examples include `'1H'`, `'13C'`, and `'E'`.
- `field` — real numeric scalar for the primary magnetic field, in tesla.
- `temperature` — non-zero real numeric scalar spin temperature, in kelvin. The input check does not require positivity.

## Calculation and outputs

The function obtains the magnetogyric ratio and multiplicity from `spin(isotope)`, creates the spin matrices with `pauli(multiplicity)`, and forms the Zeeman Hamiltonian `H=-mg_ratio*field*S.z`. Energies are returned as `E=ħ diag(H)/(k_B temperature)`, i.e. fractions of `k_B T`. The source uses exact SI constants `ħ=6.62607015e-34/(2π) J·s` and `k_B=1.380649e-23 J/K`.

To avoid overflow in the Boltzmann factors, it computes `exp(-E+min(E))` and normalises the result to obtain `P`. It returns `dP=-diff(P)`; these are signed differences for adjacent entries in the returned vector order, not absolute differences.

- `E` — vector of energy levels in units of `k_B T`.
- `P` — normalised vector of level populations.
- `dP` — signed adjacent-level population differences.

The sign of the magnetogyric ratio matters: the source notes it is negative for electrons and positive for protons. Field and temperature are checked as real scalars, with temperature additionally checked to be non-zero; the source does not check either for finiteness.

## Source link

[Spinach Wiki: levelpop.m](https://spindynamics.org/wiki/index.php?title=levelpop.m)
