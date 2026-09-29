# kernel/optimcon/distortions/szf.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/szf.m)

## Purpose and syntax

`[w,J]=szf(w,z)` applies a discrete single-zero filter to an optimal-control waveform. Each adjacent odd/even row pair stores the real and imaginary components of one complex control signal, with one time slice per column. The first sample is unchanged; for later samples the filter is `Y(k)=(X(k)-z*X(k-1))/(1-z)`.

## Inputs and output

- `w`: Real numeric waveform with an even number of rows. The output has the same dimensions as the input.
- `z`: Numeric vector with exactly one finite value per X,Y pair. Complex values are accepted; values equal to 1 are rejected, and the implementation imposes no additional magnitude bound. There is no default or scalar broadcast across multiple pairs. As a filter coefficient, `z` is dimensionless.

The source header gives a pole parametrisation `p=exp(-r*dt+1i*(omega-omega_rf)*dt)`, with damping rate `r`, pole frequency `omega`, rotating-frame frequency `omega_rf`, and time step `dt`. The executable filter consumes the supplied `z` directly; it does not compute `p` or state a formula connecting `p` to `z`.

The user is responsible for leaving sufficient ring-down margin.

## Jacobian

When requested, `J` is the Jacobian of the vectorised output with respect to the vectorised input. The implementation obtains it by automatic differentiation and returns the extracted real Jacobian.

<https://spindynamics.org/wiki/index.php?title=szf.m>
