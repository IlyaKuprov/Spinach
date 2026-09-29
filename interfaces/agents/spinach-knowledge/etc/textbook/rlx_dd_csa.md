# etc/textbook/rlx_dd_csa.m

- MATLAB implementation: [etc/textbook/rlx_dd_csa.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_dd_csa.m)

**Signature:** `[A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)`

## Purpose

Computes Redfield relaxation and cross-relaxation for two spin-1/2 particles with dipole–dipole (DD) and CSA interactions. It returns the separated CSA/DD terms as well as totals, TROSY component rates, and the longitudinal cross-relaxation rate.

## Inputs

- `B0` — finite real scalar magnetic field in tesla.
- `tau_c` — finite positive real scalar rotational correlation time in seconds.
- `isotopes` — two-element cell array of character-array isotope labels, e.g. `{'13C','19F'}`; both spins must have multiplicity 2 (spin-1/2).
- `deltas` — two-element cell array of finite, real, symmetric 3-by-3 chemical-shift tensors in ppm.
- `coords` — two-element cell array of finite real 1-by-3 Cartesian coordinate rows in ångströms.

## Calculation and outputs

For each spin, the carrier-frequency expression includes the isotropic shift, `B0*(1+trace(1e-6*delta)/3)*spin(isotope)`. The CSA tensors are scaled by `1e-6` and made traceless before their second-rank invariants are evaluated; the DD tensor is built from the two coordinates with `xyz2dd`. The isotropic-tumbling spectral-density factor is `J(omega)=tau_c/(1+omega^2*tau_c^2)`; rates combine its values at zero, each carrier frequency, and sum/difference frequencies. CSA/DD cross-invariants from `blprod` enter the TROSY cross-correlation terms.

- `A` and `B`: output structures for the first and second input spins. Each has `.r1` and `.r2` substructures containing `.csa`, `.dd`, and `.total`. Each also has `.trosy.dd`, `.trosy.csa`, and `.trosy.xc`, plus `.trosy.total_bro` and `.trosy.total_nar` for the broad and narrow TROSY doublet components. The broad/narrow values use the total DD and CSA contributions with respectively plus/minus the absolute cross-correlation contribution.
- `X`: longitudinal cross-relaxation rate between the two spins.

The source documentation does not state output-rate units; this page therefore does not assign them.

## Reference

[Spinach Wiki: rlx_dd_csa.m](https://spindynamics.org/wiki/index.php?title=rlx_dd_csa.m)
