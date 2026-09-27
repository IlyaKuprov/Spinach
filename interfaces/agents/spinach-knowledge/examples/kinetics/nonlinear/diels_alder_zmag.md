# examples/kinetics/nonlinear/diels_alder_zmag.m

- Signature: `diels_alder_zmag()`

## Purpose

Time-domain Z-magnetisation dynamics in the Diels–Alder cycloaddition of acetylene to butadiene, demonstrating the nonlinear kinetics module. Calculation time: minutes.

## Physical / mathematical content

The source imports acetylene (A), butadiene (B), and cyclohexadiene (C) spin systems from DFT outputs and adds natural-abundance ethanol (D) as a solvent subsystem. A + B → C is a second-order reaction with rate constant 25 mol/(L·s); the solvent is a spectator. Concentrations weight the initial spin state, and reaction generators couple the changing reactant/product populations to spin evolution.

## Numerical / algorithmic content

The concentration trajectory runs for 10 s in 100 LG4 steps. Interpolated concentrations set the left- and right-edge reaction generators in each two-point Lie-quadrature spin step. The script plots the Z expectation values for acetylene, butadiene, and cyclohexadiene against time; it does not perform the pulse-acquire/GPU workflow of the companion spectrum example.

## Implementation structure

- Merges the three DFT-derived species with six-spin ethanol; sets B₀ = 14.1 T and unit initial concentrations in the kinetic model.
- Initializes [A, B, C, D] to [0.01, 0.02, 0, 0.1] mol/L and excludes ethanol from the plotted species and spin observables.
- Uses 100 steps over 10 s and rate constant 25 mol/(L·s).
- Plots the three concentration-weighted species Z-magnetisation trajectories.
