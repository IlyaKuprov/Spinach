# experiments/nmr_liquids/hoesy.m

**Canonical MATLAB source:** [experiments/nmr_liquids/hoesy.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/hoesy.m) · [Spin Dynamics Wiki: hoesy.m](https://spindynamics.org/wiki/index.php?title=hoesy.m)

## Purpose

Phase-sensitive heteronuclear NOESY model. The source explicitly describes this as an ideal model; gradient and diffusion attenuation, finite-pulse losses, and experimental normalisation are outside this function.

## Input contract

Call `fid=hoesy(spin_system,parameters,H,R,K)`. `parameters.sweep` gives the two F1/F2 sweep widths in Hz, `parameters.npoints` gives the point counts for the two dimensions, and `parameters.spins` identifies the nuclei (source example: `{'15N','13C'}`). `parameters.decouple_f1` lists F1-decoupled nuclei (example: `{'1H','13C'}`); `parameters.tmix` is the mixing time in seconds; `parameters.rho0` is the initial state; and `parameters.needs` should be `{'rho_eq'}` for the thermal-equilibrium state. The source requires `sphten-liouv` and numeric, same-sized matrix inputs `H`, `R`, and `K`.

## Sequence and acquisition

Starting from `rho0`, a +90° F1 x-pulse precedes the first half of sampled F1 evolution. Each listed `decouple_f1` nucleus receives an F1 midpoint π x-pulse before the second half. The following F1 pulse creates four phase branches using ±90° x- and y-pulses; the branches are homospoiled, receive +90° F2 y-pulses, and evolve for the supplied `tmix`. The source subtracts the negative-pulse branch from the positive-pulse branch separately for the x- and y-derived components to eliminate axial peaks in F2. It then decouples the first spin entry during F2 and separately acquires cosine and sine F2 evolutions with the F2 `L+` detection state. No scalar-coupling transfer delay is specified by this function; its mixing interval is `tmix`.

## Output

`fid.cos` and `fid.sin` are the separate phase-sensitive components. Each is a two-dimensional FID with `npoints(1)` F1 evolution samples and `npoints(2)` F2 detection samples.

## References

- [HOESY sequence reference](https://doi.org/10.1021/ja00353a071)
- [HOESY sequence reference](https://doi.org/10.1039/C8CP00911B)
