# interfaces/gaussian/karplus_fit.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gaussian/karplus_fit.m) · [Wiki page](https://spindynamics.org/wiki/index.php?title=Karplus_fit.m)

**Call:** `[A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)`.

`dir_path` names a directory whose `*.log` files are Gaussian output for a dihedral scan. `atoms` is a cell array; each entry is used as a four-atom index vector. For each log, `gparse` supplies `std_geom` coordinates and a `j_couplings` matrix. The dihedral is computed from the four coordinate rows, and the corresponding coupling is `j_couplings(atoms{k}(1),atoms{k}(4))`. Thus the fitted observations are angle in degrees and scalar coupling in hertz.

The logs are parsed with `parfor`. If extraction throws for any requested atom set in a log, the catch skips that log's values for all requested sets. NaN angle/coupling entries are removed. The remaining observations are fitted by linear least squares to `A*cosd(phi)^2+B*cosd(phi)+C`, after mapping angles with `mod(phi,360)`. The function returns the three coefficients and estimated coefficient standard deviations `sA`, `sB`, and `sC`. It also plots the data and the fitted curve; it does not return fit diagnostics or the figure handle.

The standard deviations are calculated from the Studentised residual scale `sqrt(sum(residuals.^2)/(N-3))` and the Jacobian of the residual vector, estimated with `jacobianest`; the covariance expression uses the inverse of `jac'*jac`. The input check requires a character `dir_path` and a cell array `atoms`, but does not check the directory contents, vector lengths, atom-index validity, or whether the scan supports a well-conditioned fit. Other called Spinach/toolbox functions include `gparse`, `dihedral`, `jacobianest`, and the plotting helpers `kfigure`, `kgrid`, `kxlabel`, and `kylabel`.
