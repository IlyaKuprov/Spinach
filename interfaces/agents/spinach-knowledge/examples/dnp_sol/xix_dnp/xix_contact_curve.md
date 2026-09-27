# examples/dnp_sol/xix_dnp/xix_contact_curve.m

- Signature: `xix_contact_curve()`

## Purpose

Tracks transfer from electron polarization (-E_z) to proton polarization (I_z) during an X-inverse-X (XiX) DNP contact. Further information: https://doi.org/10.1021/jacs.1c09900. Calculation time: seconds.

## Physical / mathematical content

The model contains a trityl electron and two protons, with anisotropic Zeeman terms, specified coordinates, and a spin temperature of 80 K. XiX irradiation drives the electron–nuclear dynamics; the detected signal is the real expectation value of proton (L_z).

## Numerical / algorithmic content

The script evaluates a contact curve with the ESR powder simulation on the `rep_2ang_1600pts_sph` grid. It uses 80 XiX blocks, each pulse lasting 48 ns, and plots 81 samples over the corresponding total contact time.

## Implementation structure

It constructs a full Zeeman–Hilbert basis, detects proton (L_z), and calls `powder` with `@xixdnp`. The irradiation uses the electron and one proton, a 17.8 MHz electron nutation frequency, and the specified offset; the plotted curve is the real-valued result.
