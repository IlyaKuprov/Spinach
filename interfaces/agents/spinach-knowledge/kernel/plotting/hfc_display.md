# kernel/plotting/hfc_display.m

- Signature: `hfc_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Adds selected atoms' hyperfine-tensor visualisations to the current 3-D molecular axes. It draws geometry through `molplot`, then overlays tensor surfaces and optional atom labels; it does not calculate spin propagation or a spectrum.

## Tensor surfaces and coordinates

- The sampled sphere uses 31 values each for `theta` and `phi` (`k=5`, `npts=2^k-1`); `ndgrid` produces 31-by-31 grids, flattened to 1-by-961 coordinate vectors and reshaped to 31-by-31 surfaces for plotting. Surface vertices are translated by the selected atom's `props.std_geom(n,:)` coordinates.
- With `options.style='ellipsoids'`, `eig(hfc,'vector')` supplies the principal values and vectors. Each sphere coordinate is scaled by its eigenvalue and `scaling`, rotated by the eigenvectors, then translated to the atom coordinate. The ellipsoid surface is grey with face alpha 0.5. Three principal-axis lines are also drawn; positive-value lines are red and negative-value lines blue.
- The default `options.style='harmonics'` converts the tensor with `mat2sphten` and forms a real radius from rank-0, rank-1 and rank-2 spherical-harmonic terms. The resulting spherical coordinates are multiplied by `scaling` and translated to the atom. Positive radii are shaded red and negative radii blue; surface face alpha is 0.25.
- Hyperfine matrix entries are used as supplied (the source describes the eigenvalues in milliTesla); `scaling` controls their display size relative to the molecular coordinates. It is not a unit conversion.

## Options and side effects

- `options.kill_iso` defaults to false; when true, subtracts `trace(hfc)/3` from each diagonal entry before drawing.
- `options.numbers` defaults to false; `options.symbols` defaults to true.
- The function calls two `light` placements, `hold('on')` and `molplot`, then adds surfaces/text to the current axes. It finishes with square, tight, equal axis scaling; hides tick labels and tick marks; selects perspective projection and orbit camera controls; and turns the box and `kgrid` on. No colormap is selected: surface colours are supplied as RGB arrays or line colours.

## Inputs and guards

`props` must contain `std_geom`, `symbols` and `hfc.full.matrix`. `atoms` is a vector of positive integer indices within the symbol list or a cell array of character strings; source examples are `{'C','H'}` and `[1 2 5]`. `scaling` must be a positive real numeric scalar. A nonempty `conmatrix` must be logical, square, and have one row/column per atom; the documented empty-matrix path lets `molplot` infer bonds with a 1.6 Å cutoff. Ellipsoid drawing errors when the eigenvector orthogonality deviation `norm(eigvecs'*eigvecs-eye(3),2)` exceeds `1e-3`, with a message recommending harmonics. The source does not validate an unmatched atom-name cell array separately.

## Links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/hfc_display.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hfc_display.m)
