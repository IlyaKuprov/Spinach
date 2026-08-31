# interfaces/orca/oparse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/orca/oparse.m`
- Signature: `props=oparse(file_name)`
- Total lines: 479

## Purpose

A parser for ORCA text output logs, versions 2.6 to 6.1. Reads the geometry and every magnetic parameter that ORCA prints in the main output file. Syntax: props=oparse(file_name)

## Physical / mathematical content

- ORCA interfaces. They recover quantum-chemistry tensors and metadata and convert them to Spinach conventions.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `find_lines()`, `sect_end()`, `read_tensor()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 65-66: Check consistency; implemented by `grumble(file_name)`.
- Lines 68-69: Read the file, preserving the leading whitespace; implemented by `file_id=fopen(file_name,'r')`.
- Lines 73-74: Content matching is done on whitespace-trimmed lines; implemented by `log_lines=strtrim(raw_lines); nlines=numel(log_lines)`.
- Lines 77-79: ORCA rules its section banners from the first column, and indents the rulers that separate the per-nucleus blocks inside a section; implemented by `ruler_mask=~cellfun(@isempty,regexp(raw_lines,'^[-=]{5,}\s*$','once'))`.
- Lines 85-92: Property sections that list nuclei run past the banners of the sub-programs that ORCA calls while filling them in, and end only where the next property section or the final report begins; implemented by `prop_idx=banner_idx(startsWith(banner_titles, {'CARTESIAN COORDINATES','ELECTRONIC G-MATRIX','ZERO-FIELD-SPLITTING', 'ELECTRIC AND MAGNETIC HYPERFINE','CHEMICAL SHIFTS',…`.
- Lines 94-95: Read the ORCA version and refuse anything that is not an ORCA log; implemented by `version_idx=find_lines(log_lines,'^Program Version\s+\d',1,nlines)`.
- Lines 102-107: Select the parser branch: ORCA 6 renamed the hyperfine and EFG blocks from "Raw" to "Total", and the shielding section from CHEMICAL SHIFTS to CHEMICAL SHIELDINGS. The spelling used by the other branch is retained as a fallback, so that patched builds and pre-release versions are still parsed; implemented by `if version_major>=6`.
- Lines 117-118: Read the last Cartesian coordinate table in the log; implemented by `sect=find(startsWith(banner_titles,'CARTESIAN COORDINATES (ANGSTROEM)'),1,'last')`.
- Lines 138-140: ORCA decorates ghost, dummy, and capped ECP centres with punctuation that the per-nucleus property headers do not repeat; implemented by `bare_symbols=regexprep(symbols,'[^A-Za-z]','')`.
- Lines 142-149: Assign atomic numbers, leaving ghost and dummy centres at zero; implemented by `periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar', 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','A…`.
- Lines 156-159: Read the charge and the spin multiplicity from the calculation settings block; later parts of the log reprint both quantities in formats that do not always refer to the system as a whole; implemented by `charge_idx=find_lines(log_lines,'^Total Charge\s+Charge\s+\.+\s+-?\d+$',1,nlines)`.
- Lines 168-169: Read the final single point energy; implemented by `energy_idx=find_lines(log_lines,'^FINAL SINGLE POINT ENERGY',1,nlines)`.
- Lines 176-177: Read the electric dipole moment; implemented by `dipole_idx=find_lines(log_lines,'^Total Dipole Moment',1,nlines)`.
- Lines 187-188: Read the last Mulliken population analysis; implemented by `sect=find(startsWith(banner_titles,'MULLIKEN ATOMIC CHARGES'),1,'last')`.
- Lines 202-203: Spin populations are only printed for open-shell systems; implemented by `if ~all(isnan(spin_pop)), props.mulliken_spin=spin_pop; end`.
- Lines 206-209: Read the last electronic g-matrix, avoiding the individual contribution sections that ORCA prints after the total; implemented by `sect=find(ismember(banner_titles,{'ELECTRONIC G-MATRIX', 'ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN'}),1,'last')`.
- Lines 226-227: Read the last zero-field splitting tensor; implemented by `sect=find(startsWith(banner_titles,'ZERO-FIELD-SPLITTING TENSOR'),1,'last')`.
- Lines 243-244: Read the last hyperfine and electric field gradient section; implemented by `sect=find(startsWith(banner_titles,'ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE'),1,'last')`.

### Control flow inferred from the code

- Line 96: conditional branch on `isempty(version_idx)`.
- Line 107: conditional branch on `version_major>=6`.
- Line 119: conditional branch on `isempty(sect)`.
- Line 125: `while` loop over `(n<=last)&&isempty(log_lines{n}), n=n+1; end`.
- Line 126: `while` loop over `n<=last`.
- Line 128: conditional branch on `isempty(atom_line)||isnan(str2double(atom_line{2})), break; end`.
- Line 151: `for` loop over `n=1:natoms`.
- Line 153: conditional branch on `~isempty(element), props.atomic_numbers(n)=element; end`.
- Line 160: conditional branch on `~isempty(charge_idx)`.
- Line 164: conditional branch on `~isempty(mult_idx)`.
- Line 170: conditional branch on `~isempty(energy_idx)`.
- Line 178: conditional branch on `~isempty(dipole_idx)`.
- Line 181: conditional branch on `numel(numbers)>=3`.
- Line 189: conditional branch on `~isempty(sect)`.

### Key state/data transformations

- Lines 69: computes `file_id` using `file_id=fopen(file_name,'r')`.
- Lines 70: computes `orca_log` using `orca_log=textscan(file_id,'%s','delimiter','\n','whitespace','')`.
- Lines 71: computes `fclose(file_id); raw_lines` using `fclose(file_id); raw_lines=orca_log{1}`.
- Lines 74: computes `log_lines` using `log_lines=strtrim(raw_lines); nlines=numel(log_lines)`.
- Lines 75: computes `props.filename` using `props.filename=file_name`.
- Lines 79: computes `ruler_mask` using `ruler_mask=~cellfun(@isempty,regexp(raw_lines,'^[-=]{5,}\s*$','once'))`.
- Lines 80: computes `title_mask` using `title_mask=false(nlines,1)`.
- Lines 81-82: computes `title_mask(2:(nlines-1))` using `title_mask(2:(nlines-1))=ruler_mask(1:(nlines-2))&ruler_mask(3:nlines)& (~cellfun(@isempty,log_lines(2:(nlines-1))))`.
- Lines 83: computes `banner_idx` using `banner_idx=find(title_mask); banner_titles=log_lines(banner_idx)`.
- Lines 88-92: computes `prop_idx` using `prop_idx=banner_idx(startsWith(banner_titles, {'CARTESIAN COORDINATES','ELECTRONIC G-MATRIX','ZERO-FIELD-SPLITTING', 'ELECTRIC AND MAGNETIC HYPERFINE','CHEMICAL SHIFTS',…`.
- Lines 95: computes `version_idx` using `version_idx=find_lines(log_lines,'^Program Version\s+\d',1,nlines)`.
- Lines 99: computes `props.orca_version` using `props.orca_version=regexp(log_lines{version_idx(1)},'\d+(\.\d+)*','match','once')`.
- Lines 100: computes `version_major` using `version_major=str2double(regexp(props.orca_version,'^\d+','match','once'))`.
- Lines 108: computes `hfc_hdr` using `hfc_hdr='^Total HFC matrix'; hfc_alt='^Raw HFC matrix'`.
- Lines 109: computes `efg_hdr` using `efg_hdr='^Total EFG matrix'; efg_alt='^Raw EFG matrix'`.
- Lines 110: computes `shield_sect` using `shield_sect='CHEMICAL SHIELDINGS'`.
- Lines 118: computes `sect` using `sect=find(startsWith(banner_titles,'CARTESIAN COORDINATES (ANGSTROEM)'),1,'last')`.
- Lines 123: computes `first` using `first=banner_idx(sect); last=sect_end(banner_idx,nlines,first)`.

### Local helper functions

- Line 437: `find_lines()` — `function idx=find_lines(log_lines,pattern,first,last)`. Returns the last line of the section that starts at the given line, that is the line just above the next banner in the supplied index
  - Representative operation: `if first>last, idx=[]; return; end`.
  - Representative operation: `hits=regexp(log_lines(first:last),pattern,'once')`.
- Line 445: `sect_end()` — `function last=sect_end(banner_idx,nlines,first)`. Returns the first 3x3 matrix printed within the specified range of lines, located as three consecutive lines that contain nothing but
  - Representative operation: `next_banner=banner_idx(banner_idx>first)`.
  - Representative operation: `if isempty(next_banner), last=nlines; else, last=next_banner(1)-1; end`.
- Line 454: `read_tensor()` — `function [tensor,found]=read_tensor(log_lines,first,last)`.
  - Representative operation: `tensor=zeros(3,3); found=false; nrows=0`.
  - Representative operation: `for n=first:min(last,first+15)`.
- Line 469: `grumble()` — `function grumble(file_name)`. I can't think that it would be terrible of me to say -and it is occasionally true -that I need physics more than friends.
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `textscan()`, `fclose()`, `strtrim()`, `cellfun()`, `regexp()`, `false()`, `title_mask()`, `ruler_mask()`, `log_lines()`, `banner_idx()`, `startsWith()`, `find_lines()`, `version_idx()`, `str2double()`.
