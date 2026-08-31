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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 75-76: Check consistency; implemented by `if ~exist('options','var'), options={}; end`.
- Lines 79-80: Read the file; implemented by `file_id=fopen(filename,'r')`.
- Lines 84-85: Deblank all lines; implemented by `for n=1:numel(g03_output), g03_output(n)=deblank(g03_output(n)); end`.
- Lines 87-88: Refuse to process files without #p option; implemented by `for n=1:numel(g03_output)`.
- Lines 95-96: Set the default error flag; implemented by `props.error=0`.
- Lines 98-99: Set the default complete flag; implemented by `props.complete=0`.
- Lines 101-102: Parse the file; implemented by `for n=1:length(g03_output)`.
- Lines 104-105: Read charge and multiplicity; implemented by `current_line=char(g03_output(n))`.
- Lines 114-115: Read the input orientation and atomic numbers; implemented by `if strcmp(g03_output(n),'Input orientation:')`.
- Lines 124-125: Read the standard orientation and atomic numbers; implemented by `if strcmp(g03_output(n),'Standard orientation:')`.
- Lines 134-135: Read the SCF energy; implemented by `current_line=char(g03_output(n))`.
- Lines 142-143: Read the total electron spin; implemented by `if length(current_line)>4 && strcmp(current_line(1:4),'<Sx>')`.
- Lines 150-151: Read isotropic hyperfine couplings; implemented by `if strcmp(g03_output(n),'Isotropic Fermi Contact Couplings')`.
- Lines 170-171: Read, renormalise, and symmetrize anisotropic hyperfine couplings; implemented by `if strcmp(g03_output(n),'Anisotropic Spin Dipole Couplings in Principal Axis System')`.
- Lines 183-186: Renormalise the dipolar part and add the isotropic part in the laboratory basis, where it stays exactly isotropic even though Gaussian's four-decimal eigenvectors are not exactly orthonormal; implemented by `dip_part=[baa(3) bbb(3) bcc(3)]/(props.multiplicity-1)`.
- Lines 198-199: Read the g-tensor (Gaussian03 for Unix); implemented by `if strcmp(g03_output(n),'g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:')`.
- Lines 218-219: Read the g-tensor (Gaussian03 for Windows); implemented by `if strcmp(g03_output(n),'g tensor (ppm):')`.
- Lines 238-241: Read chemical shielding tensors; implemented by `if strcmp(g03_output(n),'SCF GIAO Magnetic shielding tensor (ppm):')|| strcmp(g03_output(n),'MP2 GIAO Magnetic shielding tensor (ppm):')|| strcmp(g03_output(n),'Magnetic…`.

### Control flow inferred from the code

- Line 76: conditional branch on `~exist('options','var'), options={}; end`.
- Line 85: `for` loop over `n=1:numel(g03_output), g03_output(n)=deblank(g03_output(n)); end`.
- Line 88: `for` loop over `n=1:numel(g03_output)`.
- Line 90: conditional branch on `(numel(current_line)>=2)&&strcmp(current_line(1:2),'# ')`.
- Line 102: `for` loop over `n=1:length(g03_output)`.
- Line 106: conditional branch on `length(current_line)>8 && strcmp(current_line(1:8),'Charge =')`.
- Line 115: conditional branch on `strcmp(g03_output(n),'Input orientation:')`.
- Line 117: `while` loop over `~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')`.
- Line 125: conditional branch on `strcmp(g03_output(n),'Standard orientation:')`.
- Line 127: `while` loop over `~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')`.
- Line 136: conditional branch on `length(current_line)>10 && strcmp(current_line(1:9),'SCF Done:')`.
- Line 143: conditional branch on `length(current_line)>4 && strcmp(current_line(1:4),'<Sx>')`.
- Line 151: conditional branch on `strcmp(g03_output(n),'Isotropic Fermi Contact Couplings')`.
- Line 154: `for` loop over `k=1:natoms`.

### Key state/data transformations

- Lines 80: computes `file_id` using `file_id=fopen(filename,'r')`.
- Lines 81: computes `g03_output` using `g03_output=textscan(file_id,'%s','delimiter','\n')`.
- Lines 82: computes `fclose(file_id); g03_output` using `fclose(file_id); g03_output=g03_output{1}`.
- Lines 89: computes `current_line` using `current_line=char(g03_output(n))`.
- Lines 96: computes `props.error` using `props.error=0`.
- Lines 99: computes `props.complete` using `props.complete=0`.
- Lines 107-108: computes `scan_data` using `scan_data=textscan(current_line,'%*s %*s %f %*s %*s %f', 'Delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 109: computes `props.charge` using `props.charge=scan_data{1}`.
- Lines 110: computes `props.multiplicity` using `props.multiplicity=scan_data{2}`.
- Lines 116: computes `k` using `k=n+5; m=1`.
- Lines 118: computes `S` using `S=eval(['[' char(g03_output(k)) ']']); atoms(m,:)=S(4:6); atomic_numbers(m)=S(2); k=k+1; m=m+1`.
- Lines 120: computes `natoms` using `natoms=m-1; props.inp_geom=atoms; props.natoms=natoms`.
- Lines 138: computes `props.method` using `props.method=char(scan_data{1}); props.method=props.method(3:(end-1)); props.energy=scan_data{2}`.
- Lines 146: computes `props.spin` using `props.spin=[scan_data{1} scan_data{2} scan_data{3}]; props.s_sq=scan_data{4}`.
- Lines 152: computes `props.hfc.iso` using `props.hfc.iso=zeros(natoms,1)`.
- Lines 153: computes `props.isotopes` using `props.isotopes=zeros(natoms,1)`.
- Lines 156-157: computes `first_pass` using `first_pass=textscan(S,'%*f %s %*f %*f %f %*f','Delimiter', ' ','MultipleDelimsAsOne',1)`.
- Lines 158: computes `props.hfc.iso(k)` using `props.hfc.iso(k)=first_pass{2}`.

### Local helper functions

- Line 406: `grumble()` — `function grumble(filename,options)`.
  - Representative operation: `if (~ischar(filename))||isempty(filename)`.
  - Representative operation: `error('filename must be a non-empty character string.')`.

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
