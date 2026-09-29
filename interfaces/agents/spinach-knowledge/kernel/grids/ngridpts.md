# kernel/grids/ngridpts.m

- Signature: `n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`

## Behaviour and units

Estimates a minimum spatial discretisation count for a gradient-driven experiment. `grad_amps` and `grad_durs` are matching row vectors of gradient amplitudes in T/m and durations in seconds; `isotope` is a character isotope label (for example, `'1H'`); `max_coh_order` is a signed real integer; and `sample_size` is a positive real scalar in metres.

The source sums segment magnitudes rather than allowing gradient-area cancellation: `G_eff=sum(abs(grad_amps.*grad_durs))`. It then uses `spin(isotope)` and computes `k_max=abs(max_coh_order*spin(isotope)*G_eff)`, followed by `n=ceil(k_max*sample_size/pi)`. `spin` returns the magnetogyric ratio in rad/(s*T), so `k_max` has units rad/m and the result is a dimensionless integer count. This is a worst-case spatial angular wavenumber rule, not a frequency-offset or time-evolution calculation. The function returns `n` as a scalar nonnegative integer minimum recommendation; its header cautions that several times this count may be needed for a chosen accuracy.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/ngridpts.m) · [spin.m units](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/spin.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=ngridpts.m)