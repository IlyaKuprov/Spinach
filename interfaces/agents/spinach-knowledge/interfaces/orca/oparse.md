# interfaces/orca/oparse.m

- Signature: `props=oparse(file_name)`

## Purpose

A parser for ORCA text output logs, versions 2.6 to 6.1. Reads the geometry and every magnetic parameter that ORCA prints in the main output file. Syntax: props=oparse(file_name)

## Physical / mathematical content

- ORCA interfaces. They recover quantum-chemistry tensors and metadata and convert them to Spinach conventions.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Parameters / inputs

- file_name -a character string with the file path

## Outputs

- props.filename -log file name
- props.orca_version -ORCA version string
- props.symbols -atomic symbols, 1 x natoms cell
- props.atomic_numbers -atomic numbers, 1 x natoms
- props.std_geom -atomic coordinates, natoms x 3, Angstrom
- props.natoms -number of atoms
- props.charge -total charge
- props.multiplicity -spin multiplicity
- props.energy -final single point energy, Hartree
- props.dip_moment -electric dipole moment, a.u.
- props.mulliken_chg -Mulliken atomic charges, natoms x 1
- props.mulliken_spin -Mulliken spin populations, natoms x 1
- props.g_tensor.raw -g-matrix as printed by ORCA
- props.g_tensor.matrix -symmetrised g-matrix
- props.g_tensor.eigvals -eigenvalues of the symmetrised g-matrix
- props.g_tensor.eigvecs -eigenvectors of the symmetrised g-matrix
- props.zfs.matrix -zero-field splitting tensor, cm^-1
- props.zfs.eigvals -ZFS tensor eigenvalues, cm^-1
- props.zfs.eigvecs -ZFS tensor eigenvectors
- props.hfc.full.matrix -hyperfine tensors, Gauss, natoms cell
- props.hfc.full.eigvals -hyperfine eigenvalues, Gauss, natoms cell
- props.hfc.full.eigvecs -hyperfine eigenvectors, natoms cell
- props.hfc.iso -isotropic hyperfine couplings, Gauss
- props.efg -EFG tensors, a.u.^-3, natoms cell
- props.nqi -quadrupolar tensors, Hz, natoms cell
- props.isotopes -isotopes used by ORCA, natoms cell
- props.cst -shielding tensors, ppm, natoms cell
- props.j_couplings -isotropic J-couplings, Hz, natoms x natoms
- props.chi_temps -susceptibility temperatures, K
- props.chi_tensors -molar magnetic susceptibility tensors,
- cm^3*K/mol, one cell per temperature
- Only the fields that ORCA has actually printed are returned; the
- caller should test for their presence with isfield.
- Notes: ORCA prints magnetic parameters only for the nuclei that were
- requested in the input, and labels each of them with the zero
- based index of the atom in the Cartesian coordinate table. All
- per-atom outputs above are therefore indexed by the position of
- the atom in props.std_geom, and are left empty for atoms whose
- parameters were not printed.
- When a log contains multiple geometries or multiple property
- sections, for example a geometry optimisation or a relaxed
- surface scan, the last one printed is returned.

## Implementation structure

- A parser for ORCA text output logs, versions 2.6 to 6.1. Reads the
- geometry and every magnetic parameter that ORCA prints in the main
- output file. Syntax:
- props=oparse(file_name)
- file_name -a character string with the file path
- props.filename -log file name
- props.orca_version -ORCA version string
- props.symbols -atomic symbols, 1 x natoms cell
- props.atomic_numbers -atomic numbers, 1 x natoms
- props.std_geom -atomic coordinates, natoms x 3, Angstrom
- props.natoms -number of atoms
- props.charge -total charge
