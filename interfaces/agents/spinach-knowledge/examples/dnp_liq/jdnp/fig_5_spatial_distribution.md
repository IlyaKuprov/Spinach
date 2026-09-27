# examples/dnp_liq/jdnp/fig_5_spatial_distribution.m

- Signature: `fig_5_spatial_distribution()`
- Reference: [Concilio et al., *Phys. Chem. Chem. Phys.* (2022)](https://doi.org/10.1039/d1cp04186j)
- Calculation time: seconds

## Purpose

Maps the proton DNP signal around a radical pair for the JDNP parameter choice discussed in the cited paper. The script evaluates proton polarisation after 20 ms at a grid of positions in two planes, illustrating the spatial variation of the effect.

## Model and setup

The example obtains the spin system from `system_specification()`, sets the field to 14.08 T and the microwave power to 2*pi*250e3 rad/s, and sets the microwave offset from the trityl and free-electron frequencies. It adjusts the electron-proton scalar coupling using the electron and proton Zeeman frequencies. Proton coordinates are scanned from -30 to 30 Å on 30-point grids.

## Calculation and output

For each point in the Z=0 and Y=0 planes, the script updates the proton coordinate, creates the Spinach system and basis, builds the ESR Hamiltonian and relaxation superoperator, and propagates the equilibrium state for 20 ms. The reported amplitude is the real proton Lz expectation value divided by its equilibrium value. Each plane is evaluated with a `parfor` inner loop and displayed as a colour map; the colour-bar label identifies the result as 1H DNP at 20 ms.
