# examples/kinetics/nonlinear/diels_alder_spec.m

- Signature: `diels_alder_spec()`

## Purpose

Repeated pulse-acquire experiment during the Diels–Alder cycloaddition of acetylene to butadiene, demonstrating the nonlinear kinetics module. Calculation time: hours; GPU use is hard-coded.

## Physical / mathematical content

DFT output supplies the spin systems for acetylene (A), butadiene (B), and cyclohexadiene (C); natural-abundance ethanol (D) is added as solvent. The second-order reaction A + B → C is coupled to the spin evolution, while ethanol is not a reactant. The rate constant is 25 mol/(L·s), and concentration-weighted product/reactant spin systems are evolved with the time-dependent reaction generators. Redfield/T1-T2 relaxation is configured, with nonzero rates assigned to the solvent spins.

## Numerical / algorithmic content

Concentrations are integrated for 10 s in 100 steps using the LG4 stepper, then interpolated to supply the reaction-dependent generators. The spin trajectory uses a two-point Lie quadrature. Nine pulse-acquire experiments (at integer-second time points 0–8 s) run in a `parfor` loop; each uses a GPU-resident evolution and 4096 acquired points at 4000 Hz. The collected FIDs are apodised and Fourier transformed with 16384-point zero filling for a ppm waterfall plot.

## Implementation structure

- Imports the three molecular spin systems from `acetylene.out`, `butadiene.out`, and `cyclohexadiene.out`; adds six-spin ethanol.
- Sets B₀ = 14.1 T, greedy parallelisation, and reaction rate constant 25 mol/(L·s).
- Starts at [0.01, 0.02, 0, 17.1] mol/L for A, B, C, and D; ethanol is excluded from the concentration-kinetics plot and reaction.
- Sets acquisition offset 2370 Hz, sweep 4000 Hz, 4096 points, and 16384-point zero filling.
