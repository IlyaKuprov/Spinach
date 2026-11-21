# kernel/cce/cce_spectrum.m

## Native role and invocation

Explicit correlation-spectrum transform; clexpn invokes it only when noise_correlation is used with parameters.spectrum.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `spectrum=cce_spectrum(times,correlation,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_spectrum.m)

## Purpose

Integrate a declared bath correlation into a finite-window spectrum.

## Behaviour

- Uses the explicit angular-frequency grid and rectangle or trapezoid weights with convention exp(+i*omega*t).
- A positive-time domain is extended only under an explicitly stationary correlation assumption; two-sided data are used as declared.
- Applies none or Hann window without window-power renormalisation. Optional real modulation on the original time mesh supplies a squared filter integral in seconds squared.
- One-sided reality and full-time evenness are checked independently for each correlation column, relative to that column's maximum magnitude. Independent unit rescaling and large companion columns cannot conceal invalid input; zero columns remain valid.

## Inputs and outputs

**Inputs**

- `times` — time samples in seconds.
- `correlation` — time-by-observable complex correlations.
- `parameters.spectrum` — frequency_grid (rad/s), convention=exp(+i*omega*t), sidedness=two-sided|one-sided, window=none|hann, time_domain=full|positive, weights=trapezoid|rectangle; positive-time extension additionally requires an explicit scalar logical `stationary=true`; missing, false, numeric, and nonscalar assertions are refused. Full-time data need no stationarity assertion.

**Outputs**

- `spectrum` — angular-frequency grid, transformed values, and declared convention/quadrature metadata.

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.
