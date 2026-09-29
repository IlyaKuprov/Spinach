# experiments/esr_dipolar/deer_analyt.m

- MATLAB implementation: [experiments/esr_dipolar/deer_analyt.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_analyt.m)

## Purpose

This helper evaluates an analytical DEER form factor for a two-spin model with dipolar and exchange coupling. It is a trace calculation, not a pulse-sequence simulation or a Spinach context callback.

## Interface and physical parameters

`deer=deer_analyt(D,J,t)` takes `D`, the dipolar coefficient multiplying `(1 - 3 cos²(theta)) Lz Sz` in the spin Hamiltonian, in angular-frequency units; `J`, the exchange coefficient multiplying `L*S` in the NMR convention (no factor of two), also in angular-frequency units; and a numeric real array `t` of non-negative times in seconds. The source validates `D` as a positive real scalar, `J` as a real scalar, and `t` as non-negative. The returned DEER form-factor array has the same dimensions as `t`, with value 1 at `t=0`.

## Mathematical content

The closed form uses Fresnel cosine and sine integrals with oscillatory phase set by `(D + J)t`; the prefactor and Fresnel arguments depend on `D t`. It is Kuprov's expression, cited in the source to DOI [10.1038/ncomms14842](https://doi.org/10.1038/ncomms14842). The implementation explicitly replaces the zero-time value to remove the formula's indeterminate limit.

## Scope not specified by the source

The helper does not define spin quantum numbers, an orientational distribution, pulse timings, acquisition settings, or a conversion between angular-frequency units and other conventions. It returns the form factor only; it does not perform powder averaging or fit experimental data.

Source: https://spindynamics.org/wiki/index.php?title=deer_analyt.m
