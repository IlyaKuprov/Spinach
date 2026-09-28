# examples/dnp_sol/beam_dnp/beam_contact_curve.m

- Signature: `beam_contact_curve()`

## Purpose

Plots the proton `I_z` expectation during a BEAM DNP contact, showing the transformation of electron `-E_z` into nuclear `I_z`. The experiment is described in [the associated Science Advances paper](https://doi.org/10.1126/sciadv.abq0536); the source estimates seconds to run.

## Model and calculation

At 0.3483 T (X-band), the model contains one electron and two protons, with a trityl g tensor, proton Zeeman estimates, coordinates, and spin temperature 80 K. It uses a full Zeeman-Hilbert basis and proton `Lz` detection. The powder-averaged BEAM sequence uses a `rep_2ang_800pts_sph` grid, 165 blocks, pulse durations 20.0 and 28.7 ns, 32 MHz electron nutation frequency, and the source's stated offsets and reference point. The plotted observable is the real proton-detection signal versus contact time.
