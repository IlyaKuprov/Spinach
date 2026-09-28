# examples/optimal_control/bloch_siegert/coote_badcop.m

- Signature: `coote_badcop()`

## Purpose

Reproduces the BADCOP-style selective-decoupling designs from [Coote et al.](https://doi.org/10.1038/s41467-018-05400-4), with Bloch-Siegert corrections enabled during design and validation. BADCOP1, BADCOP2, and BADCOP3 are designed and evaluated. Calculation time: minutes.

## Physical / mathematical content

- The model is a single (^{13}mathrm{C}) spin at 18.8 T, with zero scalar offset and a 53.2 ppm carrier. It builds correlated state-to-state targets over C-alpha and CO offset grids, C-beta inversion bands, and, for BADCOP2/3 only, C-beta preservation points outside the inversion band. The stated duration, RF ceilings, contraction factor, inversion bands, and carrier are taken from Table 1, the paper text, and Supplementary Figure 5.
- The paper parameters are `alpha_scale=0.91` and `pulse_dur=1e-3` s. The three designs use RF ceilings of 5.94, 4.87, and 7.22 kHz and C-beta inversion bands of 5–37, 28–35, and 10–45 ppm, respectively. The latter two also target preservation outside their inversion bands.

## Numerical / algorithmic content

- The C-alpha grid is 40–72 ppm (100 points), the CO grid 165–185 ppm (30 points), each C-beta inversion grid has 60 points, and preservation points are selected from an 80-point 5–80 ppm grid. The 1 ms pulse has 200 slices. L-BFGS runs for at most 200 iterations with the `rho_ens` ensemble correlation and BSS enabled.
- Adapted and unadapted pulses are separately optimized and then propagated with BSS physics on a 251-point 0–200 ppm validation grid; the plotted profiles compare final (M_Z), with the inversion-band boundaries marked.

## Implementation structure

- Set the single-spin model and sphten-liouv basis; construct operators, normalized states, and drift Hamiltonian; define shared paper parameters and the three variants; assemble offset-correlated targets and optimize BSS-aware pulses; optimize comparison pulses without BSS correction; validate both on the dense offset grid and plot the magnetization profiles.
