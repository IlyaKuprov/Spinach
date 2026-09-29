# kernel/plotting/efg_display.m

- Signature: `efg_display(props,atoms,scaling,conmatrix,options)`

Plots selected nuclear EFG or NQI tensors on the molecular geometry in the current figure. It returns no MATLAB output arguments.

## Inputs and selection

- `props.std_geom` and `props.symbols` are required. The geometry rows identify the atoms used by the plot; `molplot` draws the molecular framework using `conmatrix`.
- `atoms` may be a numeric array of positive integer atom indices (flattened by `atoms(:)`) or a cell array of element-symbol strings. For example, `{'N','O'}` selects atoms with those symbols and `[1 2 5]` selects those indices. A symbol selection includes every matching atom.
- For each selected atom, a nonempty `props.nqi{n}` tensor takes precedence; otherwise the routine uses nonempty `props.efg{n}`. It errors if neither is available.
- `scaling` is a positive real number. A nonempty `conmatrix` must be square with one row and column per atom.

## Styles and options

The default style is `harmonics`; the alternative is `ellipsoids`. Defaults are `options.kill_iso=false`, `options.numbers=false`, and `options.symbols=true`.

- With `kill_iso=true`, the isotropic part is removed as `efg - eye(3)*trace(efg)/3` before plotting.
- In `ellipsoids` style, the tensor eigensystem scales a sampled unit sphere along its three eigenvectors, translates it to the selected atom position, and draws it as a half-transparent grey surface. Principal-axis lines extend through the nucleus: positive eigenvalues are red and negative eigenvalues blue. This style requires an orthogonal eigensystem; the code switches to an error if `norm(V'*V-eye(3),2) > 1e-3` and recommends `harmonics`.
- In `harmonics` style, `mat2sphten` supplies ranks 0, 1, and 2, which are combined with spherical harmonics to give a real radial value `R` on the sampled sphere. Coordinates are `scaling*R*[X;Y;Z]`, translated to the atom. The signed radial value sets RGB surface colour directly: positive is half-intensity red, negative is half-intensity blue, and zero has zero RGB; the surface alpha is 0.25. No additional normalisation of `R` is applied.

The molecular framework is drawn before the tensor surfaces, and optional atom numbers and symbols are added at geometry positions. The helper also installs two lights for surface rendering. All graphical objects are added to the current figure.

## Existing syntax

`efg_display(props,atoms,scaling,conmatrix,options)`

## References

- [Source: `kernel/plotting/efg_display.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/efg_display.m)
- [Spinach Wiki: `efg_display.m`](https://spindynamics.org/wiki/index.php?title=efg_display.m)
