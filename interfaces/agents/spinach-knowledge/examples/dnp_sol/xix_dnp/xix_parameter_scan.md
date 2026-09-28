# examples/dnp_sol/xix_dnp/xix_parameter_scan.m

- Signature: `xix_parameter_scan()`

## Purpose

Maps the final proton (I_z) signal in a XiX DNP experiment over electron nutation frequency and microwave resonance offset. Further information: https://doi.org/10.1021/jacs.1c09900. Calculation time: minutes for the powder scan.

## Physical / mathematical content

The model uses a trityl electron and two protons with anisotropic Zeeman interactions, specified coordinates, and spin temperature 80 K. It reports the real final point of the proton-polarization contact curve for each parameter pair.

## Numerical / algorithmic content

The scan covers 30 nutation frequencies from 10 to 50 MHz and 120 offsets from -100 to +100 MHz, with a -13 MHz reference added to the simulated offset. Each point uses 150 XiX blocks of 48 ns pulses and the `rep_2ang_400pts_sph` powder grid. Nutation frequencies are evaluated in a `parfor` loop for each offset.

## Implementation structure

It builds a full Zeeman–Hilbert basis, detects proton (L_z), evaluates `powder`/`@xixdnp` over the two grids, and updates a contour plot with electron nutation frequency and microwave offset as its axes.
