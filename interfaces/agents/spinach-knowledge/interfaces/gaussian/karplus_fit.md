# interfaces/gaussian/karplus_fit.m

- Signature: `[A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)`

Fits Gaussian J-coupling scan data to `A*cosd(phi)^2+B*cosd(phi)+C`, where `phi` is the dihedral angle in degrees.

- `dir_path`: directory of Gaussian log files differing in the scanned dihedral.
- `atoms`: cell array of four-atom index vectors defining the dihedrals; the coupling is taken between the first and fourth atoms.
- `A,B,C`: fitted coefficients; `sA,sB,sC`: estimated standard deviations.

The function parses the logs, skips uninterpretable files, fits by linear least squares, and plots the observations and fitted curve.

[Karplus_fit.m](https://spindynamics.org/wiki/index.php?title=Karplus_fit.m)
