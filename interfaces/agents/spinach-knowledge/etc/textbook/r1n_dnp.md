# etc/textbook/r1n_dnp.m

- MATLAB implementation: [etc/textbook/r1n_dnp.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1n_dnp.m)

- Signature: `R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)`

## Purpose

Estimates nuclear longitudinal relaxation in the simple cryogenic-DNP model in the source: a bulk nuclear rate plus an electron-mediated contribution.

## Inputs

All seven arguments are required; there are no defaults.

- B0 — finite, nonzero real numeric scalar main magnetic field, in tesla.
- T — positive real numeric scalar absolute temperature, in kelvin.
- g — real numeric scalar electron g-factor; a dimensionless input, multiplied by the Bohr magneton constant in the calculation.
- T1e — positive real numeric scalar electron longitudinal relaxation time, in seconds.
- T1n_bulk — positive real numeric scalar bulk nuclear longitudinal relaxation time in seconds; `1/T1n_bulk` contributes to the relaxation rate in Hz.
- r — positive real numeric scalar electron-nuclear distance, in angstroms.
- bet — real numeric scalar angle between the magnetic field and electron-nuclear direction, in radians. The check does not restrict its range.

## Model and output

Using mu0 = 4*pi*1e-7, muB = 9.274010e-24, and kB = 1.380649e-23, the source defines

~~~matlab
sech_sq=sech(g*muB*B0/(2*kB*T))^2;
geom_dd=(1-3*cos(bet)^2)/(r/1e10)^3;
R1n=(((mu0/(4*pi))*(g*muB/B0)*geom_dd)^2)*sech_sq/T1e+1/T1n_bulk;
~~~

Here the distance is converted from angstroms to metres before the inverse-cube geometric factor is formed. R1n is documented as the nuclear relaxation rate in Hz. This is a direct evaluation of the stated simple model; the source gives no temperature range or literature citation, and its reference line is a placeholder.

## Source

[Spinach Wiki: r1n_dnp.m](https://spindynamics.org/wiki/index.php?title=r1n_dnp.m).
