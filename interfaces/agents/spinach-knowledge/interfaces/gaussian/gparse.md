# interfaces/gaussian/gparse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gaussian/gparse.m`
- Signature: `props=gparse(filename,options)`
- Total lines: 424

## Purpose

A parser for Gaussian (03, 09, 16) calculation logs. Ex- tracts all potentially useful information. Syntax: props=gparse(filename,options)

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- file_name -a character strong with a file name
- options -symmetrisation of the interaction
- tensors. By default all tensors are
- symmetrised. The symmetrisation may
- be turned off by adding the following
- strings to the options cell array:
- 'g_nosymm', 'cst_nosymm',
- 'hfc_nosymm'

## Outputs

- props.inp_geom -input geometry (Angstrom)
- props.std_geom -standard geometry (Angstrom)
- props.natoms -number of atoms
- props.method -energy method
- props.energy -SCF energy (Hartree)
- props.hfc.iso -isotropic hyperfines (Gauss)
- props.hfc.full.eigvals -HFC eigenvalues (Gauss)
- props.hfc.full.eigvecs -HFC eigenvectors
- props.hfc.full.matrix -HFC tensors (Gauss)
- props.g_tensor.eigvecs -g-tensor eigenvectors
- props.g_tensor.eigvals -g-tensor eigenvalues
- props.g_tensor.matrix -g-tensor
- props.cst -absolute shielding tensors
- props.k_couplings -isotropic K-couplings (Hz)
- props.j_couplings -isotropic J-couplings (Hz)
- props.srt -spin-rotation tensor
- props.nqi -nuclear quadrupolar tensors
- props.chi -susceptibility tensor
- props.gibbs -Gibbs free energy (Hartree)
- props.symbols -atomic symbols
- props.isotopes -nuclear isotopes used by Gaussian
- props.atomic_numbers -atomic numbers
- props.charge -overall charge
- props.el_dip_std -electric dipole moment, Debye
- props.multiplicity -overall multiplicity
- props.filename -log file name
- props.error -true if the calculation
- contains an error of any type
- Notes: the following keywords must be added to the route
- section of the Gaussian input file to produce a
- useful log:
- #p nmr=(giao,spinspin,susceptibility)
- output=pickett pop=minimal IOp(6/82=1)
- Gaussian divides its isotropic Fermi contact couplings
- by 2S=multiplicity-1, but prints the anisotropic spin
- dipole couplings without that normalisation; the two
- blocks therefore disagree by 2S for anything above a
- doublet. This is corrected here, and the hyperfine
- tensors returned are the ones that enter the spin
- Hamiltonian as S*A*I, in agreement with oparse.m

## Implementation structure

- A parser for Gaussian (03, 09, 16) calculation logs. Ex-
- tracts all potentially useful information. Syntax:
- props=gparse(filename,options)
- file_name -a character strong with a file name
- options -symmetrisation of the interaction
- tensors. By default all tensors are
- symmetrised. The symmetrisation may
- be turned off by adding the following
- strings to the options cell array:
- 'g_nosymm', 'cst_nosymm',
- 'hfc_nosymm'
- props.inp_geom -input geometry (Angstrom)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `fopen()`, `textscan()`, `fclose()`, `g03_output()`, `deblank()`, `char()`, `strcmp()`, `current_line()`, `eval()`, `atoms()`, `atomic_numbers()`, `baa()`, `bbb()`, `bcc()`.
