# interfaces/castep/c2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/castep/c2spinach.m`
- Signature: `props=c2spinach(file_name)`
- Total lines: 181

## Purpose

Parser for .magres files written by CASTEP and other codes in the CCP-NC magres v1.0 format. Reads the [atoms] and [magres] blocks and returns the geometry and the magnetic resonance tensors, keyed to the atoms in the order in which the atom records appear in the file. Syntax: props=c2spinach(file_name)

## Physical / mathematical content

- CASTEP interface. Recovers the periodic DFT (GIPAW) shielding, electric field gradient, and reduced spin-spin coupling tensors and converts them to the conventions used by `gparse` and `oparse`, so that `g2spinach` consumes the output the same way.
- Shielding tensors are returned in the printed component order (xx xy xz yx yy yz zx zy zz as the rows of the 3x3 matrix), which is the orientation in which Spinach contracts a Zeeman tensor with the spin operator on the left and the field on the right; the antisymmetric part of the shielding therefore enters with the sign the magres standard defines.
- Reduced couplings K (magres units 10^19 T^2 J^-1) are multiplied by mu_N^2/h, which is the unit `gparse` reports Gaussian K-couplings in; `g2spinach` then converts them into J-couplings for the isotopes it is given.

## Numerical / algorithmic content

- Every record is matched with a regular expression that spells out the field count and the number syntax; a record with the right tag and the wrong shape is an error rather than a silently shifted tensor.
- Block tags (`[atoms]`, `[magres]`, also the angle-bracket form) are located first; nested, unbalanced, missing, or repeated [atoms] and [magres] blocks are errors. Everything outside those two blocks, including the `[magres_old]` block CASTEP appends, is ignored.
- Units records for atom, ms, efg, and isc are checked against the standard units and anything else is an error.
- The file contains an explicit `grumble(...)` validator for the input argument.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(file_name)`.
- Lines 61-65: Read the file, drop the comments, and trim the lines; implemented by `magres_log=strtrim(regexprep(magres_log{1},'#.*$',''))`.
- Lines 67-86: Locate the block tags and refuse nested or unbalanced blocks; implemented by `tag_tokens=regexp(magres_log,'^[\[<](/?)([A-Za-z_]\w*)[\]>]$','tokens','once')` followed by a pass over the tag lines that pairs each opening tag with its closing tag and records the block name and line range; exactly one [atoms] and one [magres] block are required.
- Lines 88-97: Check the units of the records that are read below; implemented by matching `^units\s+(\S+)\s+(\S+)$` inside the two blocks against Angstrom, ppm, au, and 10^19.T^2.J^-1.
- Lines 99-100: Regular expressions for a number and for a 3x3 tensor; implemented by `num_pat='([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)'` and `ten_pat=repmat(['\s+' num_pat],1,9)`.
- Lines 102-110: Atom records: species, label, index in label, and coordinates; implemented by matching `^atom\s+(\S+)\s+(\S+)\s+(\d+)` plus three numbers; a malformed atom record or an empty atom list is an error.
- Lines 112-115: Atom keys, delimited and glued, for matching the tensor records; implemented by `atom_keys=strcat(label,{' '},index)` and `glued_keys=strcat(label,index)`; duplicate keys are an error.
- Lines 117-142: Shielding and EFG records, matched to the atoms by label and index; a record that has label, index, and nine numbers is keyed directly, a record with the label and index glued into one token (the CASTEP printing bug for indices with three or more digits) is keyed through `glued_keys`; a record that matches no atom or several atoms, or a repeated record for the same atom, is an error. `props.cst` and `props.efg` are returned only when at least one tensor was found.
- Lines 146-165: Reduced spin-spin coupling records, isotropic parts in gparse units; self-coupling records are skipped, the isotropic parts of K(A,B) and K(B,A) are averaged when both are present, and `props.k_couplings` is returned only when isc records exist.

### Local helper functions

- Line 170: `grumble()` — `function grumble(file_name)`; checks that the file name is a character string.

## Parameters / inputs

- file_name - the name of the *.magres file, a character string

## Outputs

- props.filename - log file name
- props.symbols - atomic symbols, 1 x natoms cell
- props.std_geom - atomic coordinates, natoms x 3, Angstrom
- props.natoms - number of atoms
- props.cst - chemical shielding tensors relative to the bare nucleus in vacuum, ppm, 1 x natoms cell, printed component order
- props.efg - EFG tensors, a.u., 1 x natoms cell
- props.k_couplings - isotropic reduced spin-spin couplings, natoms x natoms, in the units used by gparse.m (magres K times mu_N^2/h, Hz)
- Only the tensors that the file contains are returned; test for their presence with isfield. An atom without a tensor gets an empty cell.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `textscan()`, `fclose()`, `strtrim()`, `regexprep()`, `regexp()`, `cellfun()`, `strcmp()`, `strcat()`, `str2double()`, `reshape()`, `unique()`, `ischar()`.
- Consumers in the repository: `examples/nmr_solids/case_studies/mathies_14n_13c/*`, `examples/nmr_solids/case_studies/mathies_carbonate/*`, `examples/visualisation/efg_silicate.m` (through `efg_display`), and `g2spinach` for `cst` and `k_couplings`.
