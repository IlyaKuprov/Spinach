# kernel/plotting/molplot.m

- Signature: `molplot(xyz,conmatrix)`
- MATLAB source: [`kernel/plotting/molplot.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/molplot.m)

## Purpose and inputs

Draws molecular sticks from Cartesian coordinates; it is a geometry renderer, not a molecular simulation. `xyz` is an N-by-3 numeric real matrix in Angstroms. `conmatrix` identifies atom pairs to draw; when it is empty, the function calls `conmat(xyz,1.6)`, using the 1.6 Angstrom cutoff documented by the source.

## Rendering and side effects

For every nonzero element returned by `find(conmatrix)`, the function appends the two endpoint coordinates and a NaN separator to three row vectors. Each vector therefore has length `3*nnz(conmatrix)`. One `plot3` call draws these segments in grey, RGB `[0.5 0.5 0.5]`, with line width 1.5. A symmetric connectivity matrix enumerates both ordered entries, so the same geometric segment is represented twice.

The function has no output arguments. It checks that `xyz` is numeric, real, and has three columns; a nonempty connectivity matrix must be square and match the number of coordinate rows. It does not explicitly set axis labels, limits, aspect ratio, view, or colormap.

[Wiki reference](https://spindynamics.org/wiki/index.php?title=molplot.m)
