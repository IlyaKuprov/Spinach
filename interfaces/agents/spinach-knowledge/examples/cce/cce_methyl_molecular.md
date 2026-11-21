# cce_methyl_molecular

`cce_methyl_molecular()` plots the **36-proton APPA** curve in Fig. 7 of [Jeschke et al. (2025)](https://doi.org/10.5194/mr-6-15-2025), for the twelve methyl groups of the Finland trityl core. It does not calculate the other simulated traces in the panel. The comparison experiment is the archived bare-electron Hahn echo of FTR 1 in deuterated o-terphenyl, normalised and scaled by 0.981 as in the authors' plotting script.

The calculation is the analytical high-field limit of projected CCE-2: the proton Zeeman term and the pair Ising term commute with the isolated pair evolution, the one-proton factors are unity, and Eq. (9) supplies a factor for each of the 630 proton pairs. Intramethyl tunnelling contributes a coherent exchange term (Eq. 8), not a stochastic jump rate. The 145-node hemisphere integration follows the original nine-ring quadrature and averages the completed cluster products.

The simulation evaluates the real, normalised echo on a 2048-point grid from zero to 40 microseconds.

The labelled proton and echo CSV files reproduce the archived DFT tensors, geometry, and experimental trace (G. Jeschke et al., [Zenodo 13850793](https://doi.org/10.5281/zenodo.13850793), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). The two tunnel splittings are specified directly in the example. No archived simulation output is needed to run the example. The original three-proton disjoint-CCE illustration is retained as a dense-reference regression fixture.
