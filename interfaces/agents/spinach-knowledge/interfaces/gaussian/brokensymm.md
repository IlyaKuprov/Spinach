# interfaces/gaussian/brokensymm.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gaussian/brokensymm.m) · [Spinach Wiki: brokensymm.m](https://spindynamics.org/wiki/index.php?title=brokensymm.m)

- Signature: `J=brokensymm(props_sing,props_trip)`

## Inputs and calculation

Pass two Gaussian property structures, ordinarily the outputs of `gparse` for the singlet and triplet calculations of the same biradical. Each must contain `energy` (SCF energy in Hartree) and `s_sq` (the computed expectation value of total spin squared). The function checks for these four fields; it does not parse Gaussian log files itself.

It applies Eq. 6 of the Yamaguchi treatment cited below:

`J=(props_trip.energy-props_sing.energy)/(props_sing.s_sq-props_trip.s_sq)`

The result is then multiplied by `6.57968974479e15` to convert Hartree to Hz. The output is a scalar estimate under the Hamiltonian convention `H=-2J*(Sa.Sb)`; sign interpretation depends on retaining that convention. The source describes the estimate as order-of-magnitude and really rough, so it should not be treated as a precision exchange coupling.

The implementation depends on the input structures' numeric fields and does not perform the Gaussian calculations or additional unit conversions beyond Hartree-to-Hz.

## References

- Yamaguchi equation, Eq. 6: [doi:10.1063/1.5144696](https://doi.org/10.1063/1.5144696).
- [Spinach Wiki: brokensymm.m](https://spindynamics.org/wiki/index.php?title=brokensymm.m)
