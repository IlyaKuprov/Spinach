# examples/dnp_sol/novel_dnp/novel_parameter_scan.m

- Signature: `novel_parameter_scan()`
- Source: [`examples/dnp_sol/novel_dnp/novel_parameter_scan.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/novel_dnp/novel_parameter_scan.m)

## Purpose

Builds a two-dimensional NOVEL DNP scan of the real final-time proton longitudinal signal against electron nutation frequency and microwave resonance offset. The source defines the contact sequence as 0.25 μs, cites [Redrouthu et al., DOI 10.1063/1.5000528](https://doi.org/10.1063/1.5000528), and estimates minutes for the calculation.

## Spin model and sequence

This is a three-spin model: one electron and two 1H nuclei, with `sys.magnet=0.34` (the source calls it an X-band magnet) and `inter.temperature=80` (no unit is attached in the file). The electron Zeeman principal values are `[2.00319 2.00319 2.00258]`; the proton entries are labelled “ppm guess” and are `[0 0 5]` and `[0 5 0]`. Euler angles are entered as `(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`. Coordinates are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]`; the source gives no coordinate unit. The basis is `zeeman-hilb` with no approximation, and the detected state is the 1H Lz state.

The NOVEL setup uses spins `{'E','1H'}`, `flippulse=1`, 250 steps at `1e-9` seconds per step, the `rep_2ang_100pts_sph` powder grid, and the `aniso_eq` prerequisite. For each nutation frequency the code sets `irr_powers` to that value and computes the 90-degree pulse duration as `1/(4*irr_powers)`.

## Scan definition and algorithm

The nutation-frequency array has 30 points from `1e6` to `30e6`; the offset array has 71 points from `-35e6` to `+35e6` Hz. A reference point of `-3.3e6` is added to every offset when forming `localpar.offset=[offset+reference_point 0]`. The horizontal axis is nevertheless generated from the unshifted offset array and displayed in MHz, while the vertical axis is the nutation-frequency array displayed in MHz.

For each of the 30 nutation-frequency values, the code runs 71 powder simulations over offsets using `parfor`. Each calls `powder(spin_system,@noveldnp,localpar,'esr')`; the real part of the final contact-curve element fills `dnp_surf(n,m)`. Thus the plotted matrix is indexed by nutation frequency first and offset second. After each frequency row is completed, a 100-level filled contour plot is refreshed; it uses flat shading and the `jet` colormap.

## Result and scope

Calling the no-argument function produces the evolving contour figure; it does not return the surface as a MATLAB output. The axes and values describe this configured three-spin model, powder grid, 30-by-71 scan, and fixed reference shift. The calculation depends on Spinach, `noveldnp`, powder averaging, and MATLAB `parfor` support.
