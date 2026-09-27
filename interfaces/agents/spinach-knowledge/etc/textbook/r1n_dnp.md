# etc/textbook/r1n_dnp.m

- Signature: `R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)`

## Purpose

Estimates the nuclear longitudinal relaxation rate in a simple cryogenic DNP model that combines bulk nuclear relaxation with relaxation induced by an unpaired electron. The source does not provide a literature citation for the model.

## Model

The function adds the bulk contribution `1/T1n_bulk` to an electron-mediated term. That term uses the squared dipolar angular factor `(1-3 cos²(bet))²`, the inverse-sixth-power distance dependence, electron longitudinal relaxation time `T1e`, and the thermal factor `sech²(g μB B0/(2 kB T))`. The angle `bet` is in radians; `r` is converted from ångströms to metres in the calculation.

## Inputs

- `B0` — non-zero real main-magnet field in tesla.
- `T` — positive real absolute temperature in kelvin.
- `g` — real electron g-factor (dimensionless).
- `T1e` — positive real electron longitudinal relaxation time in seconds.
- `T1n_bulk` — positive real bulk nuclear longitudinal relaxation time in seconds.
- `r` — positive real electron–nuclear separation in ångströms.
- `bet` — real angle between the field and electron–nuclear direction, in radians.

## Output

- `R1n` — modelled nuclear longitudinal relaxation rate in Hz.

## Reference

See the [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=r1n_dnp.m). The MATLAB source contains a placeholder rather than a literature reference.
