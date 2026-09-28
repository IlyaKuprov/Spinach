# examples/dnp_sol/beam_dnp/beam_parameter_scan.m

- Signature: `beam_parameter_scan()`

## Purpose

Maps the final proton `I_z` signal from a BEAM DNP contact over microwave resonance offset and electron nutation frequency. The associated experiment is described in [Science Advances](https://doi.org/10.1126/sciadv.abq0536); the source estimates minutes for the large powder-grid calculation.

## Model and calculation

The X-band model (0.3483 T) contains an electron and two protons, with trityl g-tensor values, proton Zeeman estimates, Cartesian coordinates, and spin temperature 80 K. The BEAM sequence uses 20.0/28.7 ns pulses, 165 blocks, proton `Lz` detection, and `rep_2ang_800pts_sph` powder averaging. The script sweeps 120 offsets from −60 to +60 MHz (adding the −3.3 MHz reference) and 30 electron nutation frequencies from 20 to 40 MHz. For each offset, a `parfor` loop runs the powder calculation over nutation frequency and stores the last contact-curve point; the two-dimensional result is plotted as contours.
