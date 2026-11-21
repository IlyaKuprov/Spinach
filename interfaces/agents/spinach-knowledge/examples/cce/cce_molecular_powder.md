# cce_molecular_powder

`cce_molecular_powder()` plots the **APPA** curve in Fig. 9 of [Jeschke et al. (2025)](https://doi.org/10.5194/mr-6-15-2025). The model includes the 48 methylene protons of OX063 but no solvent protons. The experimental overlay is the archived bare-electron Hahn echo in deuterated water/glycerol with the authors' normalisation. CCE-3 and pCCE(4,8), the other simulated curves in the panel, are not included.

Under the high-field approximation, the 1128 analytical pair factors from Eq. (9) are multiplied at each of the 145 nodes of the original nine-ring hemisphere quadrature. Averaging the resulting products preserves the orientation correlations across clusters; averaging each factor first would not give the paper's powder echo.

The simulation evaluates the real, normalised echo on a 2048-point grid from zero to 40 microseconds.

The labelled proton and echo CSV files reproduce the archived DFT tensors, geometry, and experimental trace (G. Jeschke et al., [Zenodo 13850793](https://doi.org/10.5281/zenodo.13850793), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). No archived simulation output is needed to run the example. The original three-node deuteron demonstration is retained as a dense-reference regression fixture.
