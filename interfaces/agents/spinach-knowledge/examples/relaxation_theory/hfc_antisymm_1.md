# examples/relaxation_theory/hfc_antisymm_1.m

- Signature: `hfc_antisymm_1()`

## Purpose

Compares textbook longitudinal (`R1`), transverse (`R2`), and cross-relaxation (`Rx`) rates with rates extracted from Spinach’s Redfield relaxation superoperator for a two-spin system with a strongly antisymmetric hyperfine tensor.

## Physical / mathematical content

The system contains a proton and an electron at 0.33 T. Their hyperfine coupling is specified by a full, nonsymmetric 3 × 3 tensor; the example uses zero equilibrium, a 10 ps correlation time, and the lab-frame relaxation representation. The textbook rates come from `rlx_hfc`; corresponding Spinach rates are evaluated as negative expectation values of the relaxation superoperator for normalized longitudinal and transverse states. A longitudinal cross term tests the transfer rate between the two spins.

## Numerical / algorithmic content

The script constructs the spin system and an `sphten-liouv` basis with no approximation, builds the Redfield relaxation superoperator, and compares each `rlx_hfc` result against the associated matrix element of that superoperator. It also prints the complete superoperator in the IST basis. The reported calculation time is seconds.

## Implementation structure

After setting the field, isotopes, and hyperfine matrix, the script selects Redfield relaxation with zero equilibrium and lab-frame retention, then constructs the basis and relaxation superoperator. It evaluates `R1` for both spins using `Lz` states, `R2` using `L+` states, and `Rx` using a pair of `Lz` states, normalizing each state before evaluating the matrix elements. Finally, it prints the full relaxation superoperator.
