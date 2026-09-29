# kernel/pulses/iserstep.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/iserstep.m) · [Spin Dynamics Wiki: iserstep.m](https://spindynamics.org/wiki/index.php?title=iserstep.m) · [Iserles et al. paper (DOI)](http://dx.doi.org/10.1088/0305-4470/39/19/S07)

Signature: `rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)`

## Purpose and inputs

Takes one Lie-equation time step for a generator that may depend on both time and the current state. The generator handle `L(t,rho)` returns an operator in rad/s for `d_rho/d_t = -i*L(t,rho)*rho`. `t` is the start time in seconds, `dt` is the step duration in seconds, and `rho_a` and `rho_b` are the input and end-of-step state vectors. `spin_system` supplies the Spinach system used by `step`.

The second argument is a three-element cell array containing the function handle, start time, and method name. Supported methods are `'PWCL'`, `'LG2'`, `'LG4'`, `'LG4A'`, `'RKMK4'`, `'RKMK-DP5'`, `'RKMK-DP8'`, and `'RKMK-RKF45'`. The source describes `LG4` as a good balance of efficiency and numerical accuracy.

## Time sampling performed by each call

Each call advances one supplied interval `dt`; it does not create a pulse waveform or run an outer time loop. `PWCL` evaluates `L` at the left time and input state, then holds that generator constant for the step. `LG2` estimates a midpoint state by a half-step, evaluates the midpoint generator, and propagates the full step with it. `LG4` and `LG4A` use Lie-group stages within the interval. The RKMK family evaluates stage generators at tableau nodes `t+c(n)*dt` and their corresponding stage states; the named DP5, DP8, and RKF45 tableaus define those nodes and weights. The source advances one supplied `dt` per call; pulse dependence can be encoded in `L(t,rho)`, but there are no separate pulse phase/amplitude fields. It contains no plotting or file-writing/export operation.

Validation requires the three-element cell, a function handle, numeric scalar start time, a supported character method, numeric `rho_a`, and numeric scalar `dt`. The source does not impose a finite/positive check on `dt`.
