# examples/relaxation_theory/quad_relaxation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/quad_relaxation_1.m`
- Signature: `quad_relaxation_1()`
- Total lines: 51

## Purpose

14N quadrupolar relaxation in glycine in liquid state. The numerical output of Spinach is compared to the analytical equation from the textbook. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 14N quadrupolar relaxation in glycine in liquid state. The
- numerical output of Spinach is compared to the analytical
- equation from the textbook.
- Calculation time: seconds
- System specification
- Spin quantum number and quadrupolar tensor
- Relaxation theory
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Textbook relaxation rate expressions
- States of interest

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `eeqq2nqi()`, `create()`, `basis()`, `relaxation()`, `rlx_nqi()`, `state()`, `num2str()`.
