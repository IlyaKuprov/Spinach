# examples/spin_chemistry/cidnp_pumping_1.m

- Signature: `cidnp_pumping_1()`
- Source: [examples/spin_chemistry/cidnp_pumping_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/cidnp_pumping_1.m)

## Purpose

Builds the relaxation-and-pumping matrix used for an active-space projection. The source header says it simulates Equation 2 of the author's paper on chemically amplified NOEs and points to [DOI 10.1016/j.jmr.2004.01.011](https://doi.org/10.1016/j.jmr.2004.01.011); this records the source's citation, not an independent verification of the paper's results. Despite its directory and filename, this function does not evolve a radical pair or calculate a recombination yield.

## Spin system and relaxation model

- The system contains one proton and one fluorine at `sys.magnet=14.1` T. Scalar shifts are zero; the source labels the fluorine chemical-shift-anisotropy input as DFT and sets its principal values to `[-47 -16 63]` with zero Euler angles. It labels the coordinates as DFT and the scalar coupling as experimental; the coupling is `50` Hz. No coordinate or anisotropy units are stated in the script comments.
- Relaxation is Redfield with secular terms retained, `tau_c=110e-12` s, equilibrium setting `IME`, and temperature value `298` as passed to Spinach. The basis uses `sphten-liouv` with approximation `none`.
- After constructing the relaxation matrix `R`, the function forms and normalises the unit, proton-Z, fluorine-Z and product-Z operators. It applies pumping terms with strengths `1.3` on proton-Z and `34.0` on fluorine-Z through `magpump`.

## Matrix output

The active-space columns are `[U Hz Fz -HzFz]`. The function prints `P'*R*P`, identified in the source as Kuprov's matrix in Equation 2. There is no initial density operator, field sweep, time evolution or recombination kinetics in this script; describing it as a radical-pair reaction model would go beyond the code.
