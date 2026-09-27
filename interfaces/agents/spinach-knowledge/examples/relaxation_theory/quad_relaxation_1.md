# examples/relaxation_theory/quad_relaxation_1.m

- Signature: `quad_relaxation_1()`

## Purpose

Calculates longitudinal and transverse relaxation rates for liquid-state glycine's `14N` nucleus using Spinach's Redfield relaxation superoperator, and compares them with the textbook expressions returned by `rlx_nqi`.

## Physical / mathematical content

The model uses a quadrupolar interaction for a spin-1 nucleus, with quadrupole coupling `1.18e6`, asymmetry `0.53`, and a correlation time of `1e-9 s`. The relaxation superoperator is evaluated in the lab frame. Longitudinal and transverse rates are obtained from the normalized `Lz` and `L+` states as `-Lz'*R*Lz` and `-Lp'*R*Lp`, respectively, and printed beside the corresponding `rlx_nqi` values.

## Numerical / algorithmic content

The script constructs a single-`14N` spin system at `14.1 T`, converts the quadrupolar coupling and asymmetry to an NQI tensor with `eeqq2nqi`, and configures Redfield relaxation with zero equilibrium state. It uses the `sphten-liouv` formalism with no basis approximation, builds the relaxation superoperator, normalizes the two states with the 2-norm, and evaluates the four rates.
