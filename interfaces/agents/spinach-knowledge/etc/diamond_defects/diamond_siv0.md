# etc/diamond_defects/diamond_siv0.m

`[sys,inter]=diamond_siv0(parameters)`

Builds a Spinach spin-system and interactions model for the SiV0 defect. The magnetic parameters are attributed to Edmonds et al., *Physical Review B* **77**, 245205 (2008) ([doi:10.1103/PhysRevB.77.245205](https://doi.org/10.1103/PhysRevB.77.245205)).

## Parameters and model

- `parameters.silicon`: `'29Si'` adds a silicon-29 nucleus with the reported hyperfine principal values 78.9, 78.9, and 76.3 MHz; `'none'` omits silicon. Any other character isotope label adds that nucleus without a hyperfine tensor specified by this routine.
- `parameters.orientation`: `'111'`, `'110'`, or `'100'`, with the corresponding crystal plane normal aligned to the field axis.
- `parameters.n_13c`: integer from 0 to 6, selecting that many reported nearest-neighbour carbon-13 sites.

The electronic isotope label is `E3`; the g principal values are 2.0035, 2.0035, and 2.0042, and the axial zero-field splitting is 1000 MHz. Each selected carbon receives an axial hyperfine tensor with principal values 30.2, 30.2, and 66.2 MHz. The six carbon sites are arranged as two on each of three dangling-bond lines. Consequently, requesting fewer than six selects a specific subset, not an average over equivalent sites: the three lines are equivalent for the `'111'` orientation, whereas for `'110'` one line has a 54.2 MHz splitting and the other two have 30.2 MHz. To model a different isotopomer or a mixture, request all six and select or weight sites in the calling script.

The routine checks that the required fields are present, that the orientation is one of the listed character values, and that `n_13c` is an integer from 0 through 6.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_siv0.m).
