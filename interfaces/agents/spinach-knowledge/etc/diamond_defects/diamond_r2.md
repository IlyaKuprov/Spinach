# etc/diamond_defects/diamond_r2.m

`[sys,inter]=diamond_r2(parameters)`

Builds a single-electron Spinach model for the R2 self-interstitial defect. The magnetic parameters are attributed to Hunt et al., *Physical Review B* **61**, 3863 (2000) ([doi:10.1103/PhysRevB.61.3863](https://doi.org/10.1103/PhysRevB.61.3863)).

## Parameters

- `parameters.d_sign`: real scalar multiplying the axial zero-field-splitting parameter, whose magnitude is 4173 MHz.
- `parameters.orientation`: crystal plane normal, specified as `'111'`, `'110'`, or `'100'`; the chosen normal is aligned with the magnetic-field (laboratory z) axis.

The electron g principal values are 2.0019, 2.0019, and 2.0021. The routine rotates the g and zero-field-splitting tensors into the selected field orientation, then returns them in Spinach's system and interaction structures. The electron isotope label is `E3`.

The input must be a structure with both fields present; the orientation must be one of the three listed character values, and `d_sign` must be a real numeric scalar.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_r2.m).
