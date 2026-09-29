# kernel/plotting/cst_display.m

- Source: [kernel/plotting/cst_display.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/cst_display.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=cst_display.m)
- Signature: `cst_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Adds chemical-shielding-tensor surfaces at selected atomic coordinates to the current molecular plot. It supports ellipsoid and spherical-harmonic styles and returns no MATLAB output arguments.

## Inputs and defaults

- `props` — molecular data with `std_geom`, `symbols`, and per-atom `cst` tensors (as supplied by `gparse`).
- `atoms` — atom indices or element-symbol strings. Examples are `[1 2 5]` and `{'C','H'}`. Numeric indices are used in the supplied order; symbol selection finds matching atoms in geometry order.
- `scaling` — positive real scalar for tensor-surface size.
- `conmatrix` — connectivity matrix passed to `molplot`; an empty vector selects its documented 1.6 Angstrom distance-cutoff behaviour.
- `options.style` — `'harmonics'` by default, or `'ellipsoids'`.
- `options.kill_iso` — defaults to false. If true, subtracts `trace(cst)/3` times the identity before drawing.
- `options.numbers` — defaults to false; `options.symbols` defaults to true. When enabled, labels are drawn for all atoms, not only those selected for tensor surfaces.

## Rendering

The function first draws the molecule with `molplot`, samples a 31-by-31 angular grid on a unit sphere, and places each tensor surface at that atom's `props.std_geom` coordinate.

In the default `'harmonics'` style, the tensor is converted with `mat2sphten`; the rank-0, rank-1, and rank-2 coefficients multiply their corresponding spherical harmonics to form a real radial value at each sampled direction. The scaled radial surface is translated to the atom. Positive radial values are red and negative values blue; the surface has 0.25 face opacity and no mesh lines.

In `'ellipsoids'` style, the tensor is diagonalised and the sphere is scaled along the three eigen-directions, rotated by the eigenvectors, and translated to the atom. The resulting ellipsoid dimensions follow the absolute eigenvalue magnitudes. This style checks eigenvector orthogonality and errors when `norm(V'*V-eye(3),2)>1e-3`; use the harmonic style for tensors that fail that check. The ellipsoid is grey with 0.5 face opacity. Its three eigen-directions are drawn through the atom: positive eigenvalues are red, negative eigenvalues blue, and line length scales with the signed eigenvalue and `scaling`.

The current axes are set to square, tight and equal scaling with perspective projection; tick labels are hidden, and the camera orbit toolbar is enabled.

## Output

No return value; the current figure is updated with molecular geometry, tensor surfaces, eigen-axes where applicable, and optional atom labels.
