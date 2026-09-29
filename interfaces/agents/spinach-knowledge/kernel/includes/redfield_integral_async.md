# kernel/includes/redfield_integral_async.m

- Signature: script include; its nested worker has signature `brw_compute_kernel(spin_system,w,job_id,upper_lim)`.

## Purpose and scope

Evaluate the asynchronous parallel integral contribution used by the Bloch-Wangsness-Redfield and Nakajima-Zwanzig relaxation theory blocks. The include consumes coupling tensors and correlation-function data prepared by its caller; it does not define which physical interactions produced those inputs. Its result is added to the caller's relaxation superoperator `R`.

The calling theory block sets `rlx_onshell` and `rlx_shift`. With `rlx_onshell=true`, the back-rotated kernel is used and gives the Redfield form at zero shift. With `false`, the resolvent kernel is used for Nakajima-Zwanzig theory. The source documents `rlx_shift` in Hz.

## Indexing, masks, and integration

For every rank `n=1:numel(Q)`, the include visits projection indices `k,m,p,q=1:(2*n+1)`. It skips a pair unless both `cheap_norm(Q{n}{k,m})>0` and `cheap_norm(Q{n}{p,q})>0`. For retained terms, `corrfun(spin_system,n,k,m,p,q)` supplies correlation weights, rates, and a state mask for each chemical species. Only nonzero exponential weights are queued. The finite upper limit is `-1.5*(1/rate)*log(1/spin_system.tols.rlx_integration)`.

Set `A=Q{n}{k,m}` and `C=Q{n}{p,q}'`. Before masking, `B` is obtained by cleaning `L0` at `spin_system.tols.rlx_integration/abs(upper_limit)`, removing terms negligible on the integration time scale. For each species, the logical mask `states{s}` zeros rows and columns outside the active state subset in `A`, `B`, `C`, and `D`. Matrix size is inherited from the Liouville-space matrices supplied by the caller; the include imposes no fixed Hilbert- or Liouville-space dimension.

The kernel matrix is `D=B-1i*(rate-rlx_shift)*speye(size(B))` on shell, and `D=-1i*(rate-rlx_shift)*speye(size(B))` off shell. The signed contribution is `-weight*A*expmint(spin_system,B,C,D,upper_limit)`, cleaned with `1e-2*spin_system.tols.rlx_zero` before assembly.

## Parallel execution and accumulation

Each retained term gets a monotonically increasing job number. The include places `A`, `B`, `C`, and `D` in the pool `ValueStore` under per-job keys, then queues `brw_compute_kernel` with `parfeval`. A worker retrieves and removes those matrices, evaluates and cleans the integral, converts its sparse result to row/column/value triples, and stores the triples under that job's result key. The client retrieves completed futures with `fetchNext`, rethrows a worker error if present, removes each retrieved result from the store, reconstructs the sparse matrix, and adds it to `R`. If no future array was created, it reports that `Q` contained no significant elements.

## References

- [Kuprov et al., Journal of Magnetic Resonance (2011), DOI: 10.1016/j.jmr.2010.12.004](http://dx.doi.org/10.1016/j.jmr.2010.12.004)
- [Auxiliary-matrix method, DOI: 10.1063/1.4928978](http://dx.doi.org/10.1063/1.4928978)
- MATLAB source: [`kernel/includes/redfield_integral_async.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/redfield_integral_async.m)
- Existing Wiki page: https://spindynamics.org/wiki/index.php?title=redfield_integral_async.m
