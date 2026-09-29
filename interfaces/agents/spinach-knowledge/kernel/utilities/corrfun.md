# kernel/utilities/corrfun.m

## Purpose

Computes the Wigner matrix element correlation function under isotropic, axial, and rhombic rotational diffusion, returning the exponential weights, decay rates, and chemical-species state maps used in Redfield relaxation theory.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/corrfun.m>

## Behaviour

- Syntax: `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`.
- The rotational diffusion model is selected per chemical species from the number of elements in `spin_system.rlx.tau_c{s}`:
  - One correlation time: isotropic model, with `D=1/(6*tau_c)`; weight `(1/(2*n+1))*krondelta(k,p)*krondelta(m,q)` and rate `-n*(n+1)*D`.
  - Two correlation times: axial model, with `D_ax=1/(6*tau_c(1))` (rotation around the main axis) and `D_eq=1/(6*tau_c(2))` (rotation perpendicular to the main axis); weight `(1/(2*n+1))*krondelta(k,p)*krondelta(m,q)` and rate `-(n*(n+1)*D_eq+((n-m+1)^2)*(D_ax-D_eq))`.
  - Three correlation times: rhombic (fully anisotropic) model, only permitted for `n=2`, with `Dxx`, `Dyy`, `Dzz` each computed as `1/(6*tau_c(i))`. Degenerate diffusion tensors (any pair of `Dxx`, `Dyy`, `Dzz` differing by less than `1e-6*mean([Dxx Dyy Dzz])`) are rejected with an error.
- For the rhombic case, five decay rates are computed as `-(4*Dxx+Dyy+Dzz)`, `-(Dxx+4*Dyy+Dzz)`, `-(Dxx+Dyy+4*Dzz)`, `-(2*Dxx+2*Dyy+2*Dzz-2*delta)`, and `-(2*Dxx+2*Dyy+2*Dzz+2*delta)`, where `delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz)`.
- Rhombic weights use coefficients `lambda_p=sqrt(2/3)*(Dxx+Dyy-2*Dzz+2*delta)/(Dxx-Dyy)` and `lambda_m=sqrt(2/3)*(Dxx+Dyy-2*Dzz-2*delta)/(Dxx-Dyy)`, combined through a fixed coefficient matrix `h`; each weight is `(1/5)*krondelta(k,p)*h(j,m)*h(j,q)`.
- Wigner function indices are sorted in descending order: `k=[1 2 3 4 5]` in the input represents `[2 1 0 -1 -2]` for `n=2`.
- Second-rank rotational correlation times (as per Spinach input) are updated automatically if other ranks are specified.
- Input validation (`grumble`) requires the `sphten-liouv` formalism, and requires `n`, `k`, `m`, `p`, `q` to be numeric, real, integer, with `n>=0` and `k,m,p,q` in `[1,2*n+1]`.

## Inputs and outputs

Inputs:
- `spin_system` — Spinach spin system object (output of `create.m`); rotational correlation times must be supplied via `spin_system.rlx.tau_c`.
- `n,k,m,p,q` — the five indices in the ensemble-averaged Wigner function product `<D{n}{k,m}(0)*D{n}{p,q}(t)'>`.

Outputs:
- `weights` — cell array (one element per chemical species) of vectors listing the weights of the exponential components of the decays.
- `rates` — cell array (one element per chemical species) of vectors listing the decay rates (negative numbers) of the exponential components.
- `states` — cell array (one element per chemical species) of logical vectors indicating which states in the basis set belong to which chemical species.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=corrfun.m>
