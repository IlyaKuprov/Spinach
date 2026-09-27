# interfaces/gaussian/gparse.m

- Signature: `props=gparse(filename,options)`

Parses Gaussian 03, 09, or 16 logs. `filename` names an existing log; optional `options` is a cell array of strings. Tensors are symmetrized by default; `'g_nosymm'`, `'cst_nosymm'`, and `'hfc_nosymm'` disable symmetrization for the corresponding tensors.

Output fields include input and standard geometries (Å), atom identities and isotope data, SCF and Gibbs energies (Hartree), hyperfine data (Gauss), g and shielding tensors, K/J couplings and spin-rotation/quadrupolar tensors (Hz), susceptibility, and electric dipole moment (Debye). In Link1 logs, the last occurrence of each quantity is retained. Spin-rotation and quadrupolar tensors are rotated from the principal-axis frame into the standard orientation; if none is printed, the input orientation is used.

A useful Gaussian log needs `#p nmr=(giao,spinspin,susceptibility)` and `output=pickett pop=minimal IOp(6/82=1)` in its route. For multiplicity above a doublet, Gaussian divides isotropic Fermi-contact couplings by `2S=multiplicity-1` but not the anisotropic spin-dipole block; the parser reconciles that scaling before returning hyperfine data.

[Source](https://spindynamics.org/wiki/index.php?title=gparse.m)
