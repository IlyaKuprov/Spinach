# kernel/optimcon/distortions/spf.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/spf.m)

## Purpose and syntax

`[w,J]=spf(w,p)` applies a discrete single-pole filter to an optimal-control waveform. The waveform has one time slice per column; each adjacent odd/even row pair stores the real and imaginary components of one complex control signal.

For each channel, the first sample is unchanged and subsequent samples follow `Y(n)=(1-p)*X(n)+p*Y(n-1)`. The source gives `p=exp(-r*dt+1i*(omega-omega_rf)*dt)`, where `r` is the damping rate, `omega` the pole frequency, `omega_rf` the rotating-frame frequency, and `dt` the time discretisation step. The source does not prescribe a numeric unit convention for these quantities.

## Inputs and output

- `w`: Real numeric waveform with an even number of rows. It is returned with the same dimensions; columns represent time slices.
- `p`: Numeric vector with exactly one finite coefficient per X,Y pair, satisfying `abs(p(k))<1`. Complex coefficients are accepted. There is no broadcast default: a scalar is valid only when there is one pair.

The user is responsible for leaving sufficient ring-down margin.

## Jacobian

When requested, `J` is the Jacobian of the vectorised output with respect to the vectorised input. The implementation obtains it by automatic differentiation and returns the extracted real Jacobian.

<https://spindynamics.org/wiki/index.php?title=spf.m>
