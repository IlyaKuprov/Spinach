# examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m

- Signature: `fig_3_time_dep_top_row()`
- Reference: [Physical Chemistry Chemical Physics, DOI: 10.1039/D1CP04186J](https://doi.org/10.1039/d1cp04186j)
- Calculation time: seconds (per source comment)

## Purpose

Provides the top-row comparison for the JDNP time traces by reducing the `system_specification()` model to the proton and one electron. The source describes this as a demonstration that the JDNP effect vanishes when the second electron is removed; the paired bottom-row script retains both electrons.

## Calculation

The script keeps the first two isotopes, Zeeman matrices, and coordinates, replaces the scalar-coupling array with a zero-filled 2-by-2 cell array, removes the listed SRFK fields, and selects Redfield relaxation. It calculates proton (L_z) trajectories at 0.034, 0.34, and 3.4 T. At each field it sets the microwave offset from the trityl–free-electron frequency difference, propagates the thermal-equilibrium state under the ESR Hamiltonian, microwave drive, offset, and relaxation superoperator for 300 ms, and normalizes the signal by its proton equilibrium expectation value. The time step is 1 ms; each trajectory is plotted in its own field panel.
