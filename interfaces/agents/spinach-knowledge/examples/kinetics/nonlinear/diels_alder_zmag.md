# examples/kinetics/nonlinear/diels_alder_zmag.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/nonlinear/diels_alder_zmag.m

- Signature: `diels_alder_zmag()`

## Model and scope

The example follows `acetylene (A) + butadiene (B) -> cyclohexadiene (C)` and computes concentration-weighted Z magnetisation alongside the chemical concentration trajectory. A, B, and C are imported from their DFT output files; each `g2spinach` call passes 31.8 as its third argument (the script does not state its unit) and sets the minimum imported J coupling to 2.0 Hz. The source labels the six-proton, coordinate-free solvent subsystem D as natural-abundance ethanol. Spins 1-2, 3-8, 9-16, and 17-22 define A-D. The initial concentrations `[0.01, 0.02, 0, 0.1] mol/L` are supplied directly in `chem.concs` and occupy the four substance unit coordinates.

The source sets B0=14.1 T and an additive reaction record with rate 25.0 L/(mol*s). Tracing all spins with `kill_spin` gives a four-pool concentration model using the same record and `kinetics`; it is stepped with LG4 over 10 s in 100 increments. Makima interpolants without extrapolation provide A, B, and C concentrations for the spin-evolution generators; D does not participate in the reaction and is omitted from both plots. The concentration plot shows A, B, and C.

## Spin observable

The initial state combines the concentration-weighted A–C longitudinal state, scaled by the 1H level populations at 300 K, with the unit-coordinate populations. `kinetics` compiles reactants [1, 2], product 3, and the original eight-spin matching table once. Prescribed interpolated populations are embedded in unit coordinates when assembling interval-edge generators; solvent is not excited. Product unit arrival is shared equally between reactants, as required by additive closure, rather than assigned only to B. A second 100-step loop advances the magnetisation with left- and right-edge reaction generators in a two-point Lie step. The plotted traces are the real overlaps with unweighted `coil_state` Lz detection vectors for spins 1, 3, and 9, labelled acetylene, butadiene, and cyclohexadiene; the axis label is “conc.-weighted exp. value, 300K.”

The output is a model trajectory, not a reported measured magnetisation. The source header estimates minutes of calculation.
