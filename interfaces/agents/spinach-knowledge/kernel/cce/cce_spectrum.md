# kernel/cce/cce_spectrum.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_spectrum.m)

## Purpose

Integrate a declared bath correlation into a finite-window spectrum.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Uses the explicit angular-frequency grid and rectangle or trapezoid weights with convention exp(+i*omega*t).
- A positive-time domain is extended only under an explicitly stationary correlation assumption; two-sided data are used as declared.
- Applies none or Hann window without window-power renormalization. Optional real modulation on the original time mesh supplies a squared filter integral in seconds squared.

- One-sided reality and full-time evenness are checked independently for each correlation column, relative to that column's maximum magnitude. Independent unit rescaling and large companion columns cannot conceal invalid input; zero columns remain valid.

## Inputs and outputs

```matlab
spectrum=cce_spectrum(times,correlation,parameters)
```

**Inputs**

- `times` — time samples in seconds.
- `correlation` — time-by-observable complex correlations.
- `parameters.spectrum` — frequency_grid (rad/s), convention=exp(+i*omega*t), sidedness=two-sided|one-sided, window=none|hann, time_domain=full|positive, weights=trapezoid|rectangle; positive-time extension additionally requires an explicit scalar logical `stationary=true`; missing, false, numeric and nonscalar assertions are refused. Full-time data need no stationarity assertion.

**Outputs**

- `spectrum` — angular-frequency grid, transformed values and declared convention/quadrature metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_spectrum.m)
- [cce_noise.m](cce_noise.md)
