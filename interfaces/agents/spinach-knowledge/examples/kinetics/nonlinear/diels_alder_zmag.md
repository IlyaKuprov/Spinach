# examples/kinetics/nonlinear/diels_alder_zmag.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/nonlinear/diels_alder_zmag.m

- Signature: `diels_alder_zmag()`

## Model and scope

The example follows `acetylene (A) + butadiene (B) -> cyclohexadiene (C)` and computes concentration-weighted Z magnetisation alongside the chemical concentration trajectory. A, B, and C are imported from their DFT output files; each `g2spinach` call passes 31.8 as its third argument (the script does not state its unit) and sets the minimum imported J coupling to 2.0 Hz. The source labels the six-proton, coordinate-free solvent subsystem D as natural-abundance ethanol. Spins 1-2, 3-8, 9-16, and 17-22 define A-D. The merged kinetic bookkeeping concentrations are `[1, 1, 1, 1]`, while the explicit initial dynamic concentrations are `[0.01, 0.02, 0, 0.1] mol/L`.

The source sets B0=14.1 T and `rrc=25.0`; its comment labels the rate value `mol/(L*s)`, but the implemented bimolecular terms with concentrations in mol/L require `L/(mol*s)`. A nonlinear four-species generator is stepped with LG4 over 10 s in 100 increments. Makima interpolants without extrapolation provide A, B, and C concentrations for the spin-evolution generators; D does not participate in the reaction and is omitted from both plots. The concentration plot shows A, B, and C.

## Spin observable

The initial state combines each of A, B, and C's Lz state weighted by its initial concentration, then scales it using the 1H level populations at 300 K. The reaction generators are built with `react_gen` from reactants [1, 2], product 3, and an explicit eight-spin matching table. A second 100-step loop advances the magnetisation with left- and right-edge reaction generators in a two-point Lie step. The plotted traces are the real overlaps with Lz detection states for spins 1, 3, and 9, labelled acetylene, butadiene, and cyclohexadiene; the axis label is “conc.-weighted exp. value, 300K.”

The output is a model trajectory, not a reported measured magnetisation. The source header estimates minutes of calculation.
