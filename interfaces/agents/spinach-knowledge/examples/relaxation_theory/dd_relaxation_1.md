# examples/relaxation_theory/dd_relaxation_1.m

- MATLAB implementation: [examples/relaxation_theory/dd_relaxation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_relaxation_1.m)

- Signature: `dd_relaxation_1()`
- Source: [examples/relaxation_theory/dd_relaxation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_relaxation_1.m)

## Purpose and model

This short example builds a Bloch-Redfield-Wangsness relaxation superoperator for a two-spin heteronuclear pair coupled by a through-space dipolar interaction. The positions are supplied as Cartesian coordinates, so Spinach constructs the dipolar coupling from the geometry. The stochastic modulation is represented by a single correlation time, 5 ns, passed to the textbook dipolar-rate calculation; the source does not independently specify a spectral-density formula or compare alternative motional models.

## System and relaxation settings

The system is one proton and one carbon-13 at a field of 14.1 T, with coordinates [0, 0, 0] and [0, 0, 1.02]. The source does not label the coordinate unit. It selects Redfield relaxation, zero equilibrium, and lab-frame retention, and uses the full `sphten-liouv` basis (`approximation='none'`). It does not set an explicit secular restriction or select particular cross-correlation terms; no such setting should be inferred from this script.

## Rates and operators

The script constructs `R=relaxation(spin_system)` and calls `rlx_dip` for reference `R1`, `R2`, and cross-relaxation (`Rx`) rates. It also evaluates superoperator matrix elements for normalised `Lz` states on each spin (`R1`), normalised `L+` states on each spin (`R2`), and a pair of normalised `Lz` states for the transfer rate, using `-rho'*R*rho` or `-rho_b'*R*rho_a`. The displayed rate units are Hz. It prints the complete relaxation superoperator in the IST basis. This is a rate/superoperator comparison, not a time-domain signal or an experimental measurement.
