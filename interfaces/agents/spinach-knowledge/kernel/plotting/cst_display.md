# kernel/plotting/cst_display.m

- Signature: `cst_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Plots chemical shielding tensors (CSTs) for selected atoms on the molecular geometry, using ellipsoids or spherical harmonics.

## Numerical / algorithmic content

The default `harmonics` style converts each tensor to irreducible spherical-tensor coefficients and evaluates the associated spherical-harmonic surface. The `ellipsoids` style diagonalizes a symmetric tensor, scales a sampled unit sphere by its eigenvalues, rotates it by the eigenvectors, and translates it to the atom. Positive and negative values are shown in red and blue, respectively. With `kill_iso=true`, the isotropic component, `trace(cst)/3`, is subtracted before plotting.

## Syntax

```matlab
cst_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- `props` — molecular structure returned by `gparse`; must contain `std_geom`, `symbols`, and the per-atom CST data in `cst`.
- `atoms` — cell array of element symbols or vector of atom indices to display (for example, `{'C','H'}` or `[1 2 5]`).
- `scaling` — positive real factor applied to the tensor surfaces and axes.
- `conmatrix` — binary connectivity matrix; an empty value uses the 1.6 Å bond-distance cutoff described by the source.
- `options.style` — `'ellipsoids'` or `'harmonics'`; default is `'harmonics'`.
- `options.kill_iso` — remove the isotropic tensor component before plotting; default is `false`.
- `options.numbers` — show atom numbers; default is `false`.
- `options.symbols` — show atom symbols; default is `true`.

## Outputs

Updates the current figure with the molecular geometry, selected shielding tensors, and the requested atom labels; no MATLAB output arguments are returned.
