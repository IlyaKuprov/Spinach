# etc/textbook/rlx_dip.m

- MATLAB implementation: [etc/textbook/rlx_dip.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_dip.m)

**Signature:** `[r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)`

## Purpose

Calculates Redfield dipolar relaxation and longitudinal cross-relaxation rates for two spins in an isotropically tumbling liquid.

## Inputs

- `B0` — real scalar magnetic field in tesla.
- `spins` — two-element cell array of character-array isotope labels, e.g. `{'1H','15N'}`.
- `dist` — positive real scalar inter-spin distance in ångströms.
- `tau_c` — positive real scalar rotational correlation time in seconds.

## Calculation and outputs

The function places the spins at `[0 0 0]` and `[0 0 dist]` to build the DD tensor with `xyz2dd`, obtains its second-rank invariant with `blinv`, and sets `r_dif_c=1/(6*tau_c)`. It calculates each spin's `s(s+1)` factor from the isotope multiplicity and forms the Zeeman frequencies as `spin(spins{i})*B0`. The Redfield expressions use the rank-2 spectral-density values `spden(2,r_dif_c,omega)` at zero, individual Zeeman, sum, and difference frequencies.

- `r1` — two longitudinal rates in Hz, ordered as the input spins.
- `r2` — two transverse rates in Hz, ordered as the input spins.
- `rx` — longitudinal cross-relaxation rate in Hz.

Unlike `rlx_dd_csa`, the validation here does not restrict the isotopes to spin-1/2; the expressions include the computed spin-square factors.

## Reference

[Spinach Wiki: rlx_dip.m](https://spindynamics.org/wiki/index.php?title=rlx_dip.m)
