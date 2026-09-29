# examples/relaxation_theory/hfc_relaxation_1.m

- MATLAB implementation: [examples/relaxation_theory/hfc_relaxation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/hfc_relaxation_1.m)

Source: [examples/relaxation_theory/hfc_relaxation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/hfc_relaxation_1.m)

## Purpose

Constructs the full Redfield relaxation superoperator for a liquid-state proton–electron pair with a point-dipole anisotropic hyperfine interaction. It compares selected entries of that superoperator with textbook dipolar-rate expressions and prints the complete matrix; the source gives an expected run time of seconds.

## Spin system and relaxation model

The spins are 1H and E at Cartesian coordinates [0, 0, 0] and [0, 0, 1.5], in a 14.1 T field. The source describes the hyperfine coupling as derived from those coordinates using the point-dipole approximation; it does not state a coordinate unit. The liquid-motion model is Redfield relaxation with a single tau_c of 10e-12 s (10 ps), zero equilibrium, and labframe retention. The basis is the complete sphten-liouv basis (approximation none).

## Rates inspected

After forming R=relaxation(spin_system), the script obtains textbook r1, r2, and cross-relaxation rx values from rlx_dip, using the field, isotopes, coordinate separation, and correlation time. It then normalises longitudinal Lz states and compares their -rho' * R * rho values with the two spins' printed R1 rates; transverse L+ states are used similarly for R2. The cross term is evaluated between normalised longitudinal states as -rho_b' * R * rho_a. Printed rate labels are Hz. Finally, the complete R is displayed in the IST basis.

This is a model-level comparison between the Redfield matrix elements and the textbook calculation in the script. The source contains no recorded rate output or experimental measurement, so this description makes no numerical agreement or experimental-validation claim.
