# kernel/includes/redfield_integral_serial.m

- Signature: script include; it has no independent function signature.

## Purpose and scope

Evaluate the serial integral contribution used by the Bloch-Wangsness-Redfield and Nakajima-Zwanzig relaxation theory blocks. It consumes caller-prepared coupling tensors and correlation data and adds each evaluated contribution to the relaxation superoperator `R`; it does not define which physical interactions produced those inputs.

The calling theory block sets `rlx_onshell` and `rlx_shift`. On shell, the back-rotated kernel gives the Redfield form at zero shift. Off shell, the resolvent kernel corresponds to Nakajima-Zwanzig theory. The source documents `rlx_shift` in Hz.

## Indexing, masks, and integration

For each rank `n=1:numel(Q)`, loop over projection indices `k,m,p,q=1:(2*n+1)`. The include proceeds only if `cheap_norm(Q{n}{k,m})>0` and `cheap_norm(Q{n}{p,q})>0`. `corrfun(spin_system,n,k,m,p,q)` supplies weights, rates, and state masks by chemical species; zero-weight exponential terms are skipped. The integration limit is `-1.5*(1/rate)*log(1/spin_system.tols.rlx_integration)`.

The coupling matrices are `A=Q{n}{k,m}` and `C=Q{n}{p,q}'`. Before masking, `B` is obtained by cleaning `L0` at `spin_system.tols.rlx_integration/abs(upper_limit)`, removing terms negligible on the integration time scale. The state mask `states{s}` zeros rows and columns outside the active subset of `A`, `B`, `C`, and `D`. Their dimensions come from the caller's Liouville-space matrices; no fixed system dimension is imposed here. The kernel is `D=B-1i*(rate-rlx_shift)*speye(size(B))` on shell and `D=-1i*(rate-rlx_shift)*speye(size(B))` off shell. Each signed term is `-weight*A*expmint(spin_system,B,C,D,upper_limit)`; it is cleaned using `1e-2*spin_system.tols.rlx_zero` and accumulated into `R`.

## Execution path

This version evaluates terms directly in the current execution context and accumulates them serially; it does not queue futures or use the pool `ValueStore`. After the loops it clears intermediate variables and reports completion of the relaxation-integral evaluation.

## References

- [Kuprov et al., Journal of Magnetic Resonance (2011), DOI: 10.1016/j.jmr.2010.12.004](http://dx.doi.org/10.1016/j.jmr.2010.12.004)
- [Auxiliary-matrix method, DOI: 10.1063/1.4928978](http://dx.doi.org/10.1063/1.4928978)
- MATLAB source: [`kernel/includes/redfield_integral_serial.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/redfield_integral_serial.m)
- Existing Wiki page: https://spindynamics.org/wiki/index.php?title=redfield_integral_serial.m
