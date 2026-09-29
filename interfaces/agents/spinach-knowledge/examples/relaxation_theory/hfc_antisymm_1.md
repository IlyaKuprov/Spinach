# examples/relaxation_theory/hfc_antisymm_1.m

- MATLAB implementation: [examples/relaxation_theory/hfc_antisymm_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/hfc_antisymm_1.m)

- Signature: `hfc_antisymm_1()`
- Source: [examples/relaxation_theory/hfc_antisymm_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/hfc_antisymm_1.m)

## Purpose and model

This example compares textbook and Redfield-superoperator longitudinal, transverse, and cross-relaxation rates for a proton-electron pair whose hyperfine tensor has a substantial antisymmetric component. The tensor is `1e6` times `[[10, 1, 1.5], [2, 0, 3], [2.5, 1, -3]]`; the script multiplies it by `2*pi` when passing it to `rlx_hfc`. The tensor's unit is not stated on that assignment. The script passes a single 10 ps correlation time to `rlx_hfc`, but does not itself state a spectral-density function. The displayed rates are in Hz.

## System and relaxation settings

The spins are `1H` and `E` at `sys.magnet=0.33`, with one correlation time of 10 ps. The source selects Redfield relaxation, zero equilibrium, lab-frame retention, and the full `sphten-liouv` basis (`approximation='none'`). It does not set an explicit secular restriction or select particular cross-correlations, so the page does not attribute either choice to the example.

## Rates and operators

After forming `R=relaxation(spin_system)`, the script calls `rlx_hfc` for textbook rates. It evaluates normalised `Lz` state matrix elements for `R1` on each spin, normalised `L+` state matrix elements for `R2` on each spin, and a normalised pair of `Lz` states for the transfer rate `Rx`; the Redfield values are formed as `-rho'*R*rho` or `-rho_b'*R*rho_a`. It prints these rate comparisons and the complete relaxation superoperator in the IST basis. These operator matrix elements are a superoperator diagnostic, not a simulated or observed experimental signal.
