# examples/nqr/nutation_nqr_iodine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nqr/nutation_nqr_iodine.m`
- Signature: `nutation_nqr_iodine()`
- Total lines: 69

## Purpose

Powder NQR nutation curve for a system with a single 127I nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder NQR nutation curve for a system with a
- single 127I nucleus.
- Calculation time: seconds
- System specification
- Formalism and basis
- Relaxation theory
- Spinach housekeeping
- Experiment parameters
- Get a figure started
- Loop over the pulse durations
- Set pulse duration
- Run the simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `operator()`, `kfigure()`, `scale_figure()`, `powder()`, `subplot()`, `ktitle()`, `num2str()`, `kxlabel()`, `ylim()`, `set()`.
