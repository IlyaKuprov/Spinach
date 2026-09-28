# kernel/plotting/efg_display.m

- Signature: `efg_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Plots electric-field-gradient (EFG) or nuclear-quadrupole-interaction (NQI) tensors for selected atoms on the molecular geometry, using ellipsoids or spherical harmonics.

## Numerical / algorithmic content

For each selected nucleus, the routine uses its NQI tensor when present, otherwise its EFG tensor, and errors if neither is available. The default `harmonics` style converts the tensor to irreducible spherical-tensor coefficients and evaluates the spherical-harmonic surface. The `ellipsoids` style diagonalizes a symmetric tensor, scales a sampled unit sphere by its eigenvalues, rotates it by the eigenvectors, and translates it to the atom. Surface sign is coloured red for positive and blue for negative values. With `kill_iso=true`, the isotropic component, `trace(efg)/3`, is subtracted before plotting.

## Syntax

```matlab
efg_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- `props` — structure from `c2spinach` or `gparse`, with molecular geometry and symbols and per-atom `nqi` or `efg` data.
- `atoms` — cell array of element symbols or vector of atom indices to display (for example, `{'N','O'}` or `[1 2 5]`).
- `scaling` — positive real factor applied to the tensor surfaces and axes.
- `conmatrix` — binary connectivity matrix; an empty value uses the 1.6 Å bond-distance cutoff described by the source.
- `options.style` — `'ellipsoids'` or `'harmonics'`; default is `'harmonics'`.
- `options.kill_iso` — remove the isotropic tensor component before plotting; default is `false`.
- `options.numbers` — show atom numbers; default is `false`.
- `options.symbols` — show atom symbols; default is `true`.

## Outputs

Updates the current figure with the molecular geometry, selected tensors, and requested atom labels; no MATLAB output arguments are returned.
