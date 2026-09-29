# interfaces/orca/oparse.m

Source: [interfaces/orca/oparse.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/orca/oparse.m)
Wiki: [oparse.m](https://spindynamics.org/wiki/index.php?title=oparse.m)

## Interface

`props = oparse(file_name)` reads one ORCA main-output text log and returns a structure. The source documents ORCA versions 2.6 through 6.1. `file_name` must be a MATLAB character array; the implementation rejects non-character inputs. It requires an ORCA version banner and a Cartesian-coordinate table. A log written with ORCA's miniprint option, without the atom-resolved coordinate table, cannot be imported.

The parser keeps the last Cartesian coordinate table and the last occurrence of each recognised property section, so optimisation or scan logs yield their final printed geometry/properties. It branches on the ORCA major version for section names changed in ORCA 6 (for example, raw versus total HFC/EFG blocks and chemical shifts versus shieldings), retaining fallback headings. It uses the atom order in the selected coordinate table as the indexing basis for per-atom fields. Magnetic properties are not guaranteed: ORCA prints them only for nuclei requested in the calculation, and atom-resolved entries not printed remain empty. Test optional fields with `isfield`.

## Returned data

- `filename`, `orca_version`: input path and parsed version string.
- `symbols` (1 x natoms cell), `atomic_numbers` (1 x natoms), `std_geom` (natoms x 3, Angstrom), and `natoms`: atoms and Cartesian geometry. Ghost/dummy centres have atomic number zero.
- `charge`, `multiplicity`, `energy`: total charge, spin multiplicity, and final single-point energy in Hartree.
- `dip_moment`: three-component electric dipole moment in atomic units.
- `mulliken_chg`, `mulliken_spin`: atomic charge and spin-population vectors (natoms x 1); the spin field is present only when spin populations are printed.
- `g_tensor.raw`, `g_tensor.matrix`, `g_tensor.eigvals`, `g_tensor.eigvecs`: the printed g matrix, its symmetrised matrix, and that matrix's eigenvalues/eigenvectors.
- `zfs.matrix`, `zfs.eigvals`, `zfs.eigvecs`: zero-field-splitting tensor data; matrix/eigenvalues are in cm^-1.
- `hfc.full.matrix`, `hfc.full.eigvals`, `hfc.full.eigvecs`: per-atom hyperfine tensors and decompositions in Gauss, stored in natoms-length cell arrays. `hfc.iso` contains isotropic hyperfine couplings in Gauss and uses NaN where not printed.
- `efg`, `nqi`, `isotopes`: per-atom cell arrays for electric-field-gradient tensors (a.u.^-3), quadrupolar tensors (Hz), and ORCA's isotope assignments.
- `cst`: per-atom shielding tensors in ppm.
- `j_couplings`: isotropic J-couplings in Hz (natoms x natoms).
- `chi_temps`, `chi_tensors`: susceptibility temperatures in K and one molar-susceptibility tensor per temperature (cm^3*K/mol).

Fields are added only when the corresponding section/data are found; the geometry and version are required. Tensor outputs follow the parser's 3 x 3 matrix decomposition, and per-nucleus tensor properties are held in atom-indexed cells.
