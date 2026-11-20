# Barros et al. (2009), Figure 2: digitised experimental points

Source: W. Barros Jr., D. F. Gochberg, and J. C. Gore, “Nuclear magnetic resonance signal dynamics of liquids in the presence of distant dipolar fields, revisited,” *J. Chem. Phys.* **130**, 174506 (2009), [DOI: 10.1063/1.3116107](https://doi.org/10.1063/1.3116107), [PMC2832049](https://pmc.ncbi.nlm.nih.gov/articles/PMC2832049/). The source was the [442 × 715 PMC figure JPEG](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/2578/2832049/4d8fb61878d5/JCPSA6-000130-174506_1-g002.jpg), SHA-256 `b93a2d0b1e67ff175b5c0011480a735dd55f3a9d6d430c458116083d0a65e71d`.

[barros_fig2.csv](../../../../../../examples/imaging/distant_dipolar/barros_fig2.csv) contains 665 resolved open-circle centres. Panel A is parallel to the polarising field and panel B is perpendicular. In both panels the black, orange/red, blue, and green experimental series correspond to 0.05, 0.10, 0.15, and 0.20 T/m, respectively (Figure 2 caption and Section 6: the quoted inverse diffusion parameter scales as the inverse square of gradient strength). Signal is the plotted amplitude as a fraction of M₀, not re-normalised separately by curve. Only symbols were measured; no points were taken from the solid Eq. (1) curves.

| Panel | 0.05 T/m | 0.10 T/m | 0.15 T/m | 0.20 T/m |
|---|---:|---:|---:|---:|
| A | 78 | 91 | 95 | 88 |
| B | 94 | 96 | 64 | 59 |

## Calibration and extraction

Least-squares affine calibration used the labelled 0–0.6 s major ticks on the top frame and 0–0.12 signal major ticks on the right frame; zero-based pixel centres are `x` to the right and `y` down. Forward calibration is `x = x0 + sx*t2_ms` and `y = y0 + sy*signal`.

| Panel | x0 (px) | sx (px/ms) | y0 (px) | sy (px per unit amplitude) | Largest tick residual, x / y (px) |
|---|---:|---:|---:|---:|---:|
| A | 54.76786 | 0.5826786 | 301.14286 | −2310.7143 | 0.30 / 0.29 |
| B | 55.50000 | 0.5914286 | 674.48214 | −2318.7500 | 0.21 / 0.39 |

Two separate pixel measurements were made in colour-specific search corridors. Where the circles are densely spaced, pass 1 locates bright-aperture extrema using a difference of Gaussian smoothings (σ = 0.48 and 1.4 px). Pass 2 independently locates minima in ink relative to both opposing tangential flanks (±1.65 px along the marker band), averaged over three short normal offsets. Local parabolic peak positions were paired one to one within 1.6 px. For the clearly open isolated circles near the start, pass 1 uses a thresholded hole centroid and pass 2 uses a radially weighted ring-versus-centre contrast peak, paired within 2 px. A source-image overlay was inspected to exclude background edges, neighbouring colours, and unresolved or occluded symbols. Apparent marker pitch was used only to reject artefacts in dense stretches. Neither time nor signal was snapped to a grid, inferred from a curve, or interpolated at missing symbols. Overlapping low-signal tails are only partially represented, particularly for the blue and green series in panel B.

The main `t2_ms,signal` columns are the arithmetic means of the corresponding `*_pass1*` and `*_pass2*` columns. The uncertainty columns are practical digitisation half-widths, not experimental uncertainties or statistical confidence intervals. Each adds a 0.5 px calibration allowance, a 1 px symbol-centre allowance, and half the measured pass separation, component-wise after axis conversion. Median absolute pass separations were 0.058 ms and 0.000047 in signal; maximum separations were 2.061 ms and 0.000500. The time half-widths range from 2.537 to 3.605 ms, and the signal half-widths from 0.000647 to 0.000899. Passes share the same JPEG and axis calibration, so their agreement cannot quantify shared systematic error.
