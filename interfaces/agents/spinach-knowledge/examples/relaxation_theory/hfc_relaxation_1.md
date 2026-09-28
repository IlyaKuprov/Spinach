# examples/relaxation_theory/hfc_relaxation_1.m

- Signature: `hfc_relaxation_1()`

## Purpose

Builds the Redfield relaxation superoperator for a liquid-state proton–electron pair whose anisotropic hyperfine coupling is obtained from the point-dipole approximation, then compares its rate matrix elements with textbook rates.

## Physical / mathematical content

The two spins are placed at `[0 0 0]` and `[0 0 1.5]`, with a 14.1 T field. Redfield relaxation uses zero equilibrium, lab-frame retention, and a 10 ps correlation time. The script obtains longitudinal, transverse, and cross-relaxation rates from `rlx_dip`, using the inter-spin distance computed from the coordinates.

## Numerical / algorithmic content

Spinach constructs an `sphten-liouv` basis without approximation and assembles the relaxation superoperator. Normalized `Lz` states provide the `R1` matrix elements, normalized `L+` states provide `R2`, and a pair of `Lz` states gives the cross-relaxation element. Each is displayed alongside the corresponding `rlx_dip` result in Hz; the full relaxation superoperator is then printed in the IST basis. The calculation time is seconds.

## Implementation structure

The script specifies field, isotopes, and Cartesian coordinates; configures the Redfield model; creates the spin system and basis; and calls `relaxation`. It computes textbook rates with `rlx_dip` and the coordinate separation, evaluates the matching Spinach rates for both spins, compares the longitudinal cross term, and prints `full(R)`.
