# kernel/plotting/hfc_display.m

- Signature: `hfc_display(props,atoms,scaling,conmatrix,options)`

## Purpose

Displays selected atoms' hyperfine tensors in a 3D molecular figure. `options.style` selects ellipsoid or spherical-harmonic surfaces; the default is `harmonics`.

## Physical / mathematical content

- For ellipsoids, the tensor's principal values set the axis lengths and its eigenvectors set their directions; plotted principal-axis colours distinguish positive and negative values.

## Numerical / algorithmic content

- The function draws the molecular geometry, selects requested atoms, and plots their hyperfine tensors using the chosen style. `options.kill_iso` can remove the isotropic tensor component before plotting.
- Tensor diagonalisation is used to orient and scale ellipsoids; this is a visualization step, not a spectral or stationary-state calculation.

## Syntax

```matlab
hfc_display(props,atoms,scaling,conmatrix,options)
```

## Parameters / inputs

- `props` - output of `gparse.m`
- `atoms` - cell array of element symbols or a vector of atom indices, e.g. `{'C','H'}` or `[1 2 5]`
- `scaling` - positive real scalar factor for scaling the tensors in the visualization
- `conmatrix` - binary connectivity matrix; if empty, a 1.6 Å distance cutoff is used
- `options.style` - `ellipsoids` or `harmonics` (default: `harmonics`)
- `options.kill_iso` - set true to remove isotropic tensor components before plotting (default: false)
- `options.numbers` - set true to display atom numbers (default: false)
- `options.symbols` - set false to hide atom symbols (default: true)
