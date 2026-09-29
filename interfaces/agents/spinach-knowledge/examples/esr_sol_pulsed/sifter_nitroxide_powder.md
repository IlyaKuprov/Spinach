# examples/esr_sol_pulsed/sifter_nitroxide_powder.m

- Signature: `sifter_nitroxide_powder()`

## Purpose

Builds a powder-averaged SIFTER example and displays the two-dimensional signal and its diagonal. The source comments estimate calculation time as minutes.

## Spin system and sequence setup

The four-spin ordering is two electrons followed by two `14N` nuclei. The electron Zeeman principal values are both `[2.0087, 2.0058, 2.0018]`, with zero Euler angles. The electron coordinate entries are `[0,0,0]` and `[0,0,20]`; the source labels them as coordinates for inter-electron DD but does not annotate their units. The nitrogen coordinate entries are empty. The coupling tensor entries assigned are for pairs (1,3) and (2,4), each with principal values `[19.8977, 20.1780, 102.8516] * 1e6` and zero Euler angles. The source does not state a unit for these coupling values. The nitrogen spins are placed in the longitudinal basis; no relaxation parameters are set in this example.

The field parameter is set to `sys.magnet=0.33`; the source gives no unit annotation. The basis is `sphten-liouv` with no approximation. The SIFTER setup starts from electron `Lz`, detects with electron `L+`, uses electron `Lx` and `Ly` pulse operators, and sets the offset to zero.

## Powder signal and figures

The calculation is `imag(powder(spin_system,@sifter,parameters,'esr'))`, with 200 points, an 8 ns timestep, and the `rep_2ang_3200pts_sph` grid. The plotted time axis is generated from the configured timestep and point count and displayed in ns. The first panel shows the 2D SIFTER signal as an image; the second plots its diagonal against time. The plotted diagonal is a view of this calculated matrix, not an independently fitted or measured observable.

## Scope

The file specifies the spin-system entries and sequence inputs for this illustrative calculation. It provides no experimental comparison or numerical reference signal, and does not assign physical units to the field, coordinate, or coupling entries; the plotted time axis is explicitly displayed in ns.

Source code: [`examples/esr_sol_pulsed/sifter_nitroxide_powder.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/sifter_nitroxide_powder.m).
