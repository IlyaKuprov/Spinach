# examples/optimal_control/bloch_siegert/coote_badcop.m

- Signature: `coote_badcop()`
- Source: [examples/optimal_control/bloch_siegert/coote_badcop.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/bloch_siegert/coote_badcop.m)
- Publication cited in the source: [Coote et al. DOI 10.1038/s41467-018-05400-4](https://doi.org/10.1038/s41467-018-05400-4)

## Objective and source-linked parameter set

The source estimates calculation time in minutes. Estimated calculation time: minutes. The example formulates BADCOP-style selective carbon decoupling pulse design with Bloch-Siegert (BSS) corrections. Its code comments associate durations, RF ceilings, contraction factor, and inversion bands with Table 1 and the article text, and the carrier with Supplementary Figure 5 of the cited paper.

The model is a single `13C` spin at 18.8 T (approximately 800 MHz for `1H`) with zero scalar offset at the reference carrier. The carrier is 53.2 ppm, the contraction factor `alpha_scale` is 0.91, and total pulse duration is 1 ms. The optimiser uses 200 equal slices (5 microseconds each), L-BFGS with at most 200 iterations, and a correlated ensemble objective. The three RF ceilings are 5.94 kHz (BADCOP1), 4.87 kHz (BADCOP2), and 7.22 kHz (BADCOP3); each is converted to an angular-frequency control level as `2*pi*rf_hz`.

## Offset targets and constraints

The offset list combines 100 C-alpha points from 40 to 72 ppm, 30 C-prime points from 165 to 185 ppm, and 60 C-beta inversion points within each variant's band: 5-37 ppm for BADCOP1, 28-35 ppm for BADCOP2, and 10-45 ppm for BADCOP3. For C-alpha offsets, alternating `Ix` and `Iy` states target their free precession under the offset Hamiltonian for `alpha_scale * pulse_duration`. C-prime and C-beta inversion points map `Iz -> -Iz`.

Only BADCOP2 and BADCOP3 add C-beta preservation targets. Their additional 80-point grid spans 5-80 ppm, excluding the corresponding inversion band; each such point maps `Iz -> Iz`. The optimiser uses the chemical-shift offsets with `Lz` as the offset operator and has BSS corrections enabled for the corrected design.

## Evaluation observable

The code also optimises a BSS-disabled waveform for each variant. It evaluates both waveforms with BSS physics present over 251 points from 0 to 200 ppm, propagates with `shaped_pulse_xy` and `expv-pwc`, and records final longitudinal magnetisation as `real(Iz' * rho)`. The plotted profiles mark the specified C-beta inversion-band edges.
