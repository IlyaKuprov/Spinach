# interfaces/castep/c2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/castep/c2spinach.m`
- Signature: `props=c2spinach(file_name)`
- Total lines: 102

## Purpose

Reads the "new format" section of CASTEP .magres files. Syntax: props=c2spinach(file_name)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(file_name)`.
- Lines 37-38: Read the file; implemented by `file_id=fopen(file_name,'r')`.
- Lines 43-44: Deblank all lines; implemented by `for n=1:numel(castep_log)`.
- Lines 48-49: Cartesian coordinates; implemented by `props.symbols={}; props.std_geom=zeros(0,3)`.
- Lines 59-60: Shielding tensors; implemented by `props.cst={}`.
- Lines 72-73: Electric field gradients; implemented by `props.efg={}`.
- Lines 85-86: Atom count; implemented by `props.natoms=size(props.std_geom,1)`.

### Control flow inferred from the code

- Line 44: `for` loop over `n=1:numel(castep_log)`.
- Line 50: `for` loop over `n=1:numel(castep_log)`.
- Line 51: conditional branch on `(numel(castep_log{n})>3)&&strcmp(castep_log{n}(1:4),'atom')`.
- Line 61: `for` loop over `n=1:numel(castep_log)`.
- Line 62: conditional branch on `(numel(castep_log{n})>1)&&strcmp(castep_log{n}(1:2),'ms')`.
- Line 65: conditional branch on `isempty(cst_spec{end})`.
- Line 74: `for` loop over `n=1:numel(castep_log)`.
- Line 75: conditional branch on `(numel(castep_log{n})>2)&&strcmp(castep_log{n}(1:3),'efg')`.
- Line 78: conditional branch on `isempty(efg_spec{end})`.

### Key state/data transformations

- Lines 38: computes `file_id` using `file_id=fopen(file_name,'r')`.
- Lines 39: computes `castep_log` using `castep_log=textscan(file_id,'%s','delimiter','\n')`.
- Lines 40: computes `fclose(file_id); castep_log` using `fclose(file_id); castep_log=castep_log{1}`.
- Lines 41: computes `props.filename` using `props.filename=file_name`.
- Lines 45: computes `castep_log(n)` using `castep_log(n)=deblank(castep_log(n))`.
- Lines 49: computes `props.symbols` using `props.symbols={}; props.std_geom=zeros(0,3)`.
- Lines 52-53: computes `atom_spec` using `atom_spec=textscan(castep_log{n},'atom %s %s %f %f %f %f', 'Delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 54: computes `props.symbols{end+1}` using `props.symbols{end+1}=atom_spec{1}{1}`.
- Lines 55: computes `props.std_geom(end+1,:)` using `props.std_geom(end+1,:)=[atom_spec{4:6}]`.
- Lines 60: computes `props.cst` using `props.cst={}`.
- Lines 63-64: computes `cst_spec` using `cst_spec=textscan(castep_log{n},'ms %s %f %f %f %f %f %f %f %f %f %f', 'Delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 68: computes `props.cst{end+1}` using `props.cst{end+1}=reshape([cst_spec{3:11}],[3 3])`.
- Lines 73: computes `props.efg` using `props.efg={}`.
- Lines 76-77: computes `efg_spec` using `efg_spec=textscan(castep_log{n},'efg %s %f %f %f %f %f %f %f %f %f %f', 'Delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 81: computes `props.efg{end+1}` using `props.efg{end+1}=reshape([efg_spec{3:11}],[3 3])`.
- Lines 86: computes `props.natoms` using `props.natoms=size(props.std_geom,1)`.

### Local helper functions

- Line 91: `grumble()` — `function grumble(file_name)`. Malevolence lurks in men who avoid wine, games, the company of beautiful women, and conversations at dinner. Such people
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

## Parameters / inputs

- file_name -the name of the *.magres file, a
- character string

## Outputs

- props.std_geom -atomic coordinates (Angstrom)
- props.symbols -atomic symbols
- props.natoms -number of atoms
- props.filename -log file name
- props.cst -chemical shielding tensors relative
- to the bare nucleus in vacuum, ppm
- props.efg -EFG tensors, a.u.^-3
- Notes: CASTEP has a printing bug with over 100 nuclei. A kind of
- workaround was implemented, but keep an eye on it.

## Implementation structure

- Reads the "new format" section of CASTEP .magres files. Syntax:
- props=c2spinach(file_name)
- file_name -the name of the *.magres file, a
- character string
- props.std_geom -atomic coordinates (Angstrom)
- props.symbols -atomic symbols
- props.natoms -number of atoms
- props.filename -log file name
- props.cst -chemical shielding tensors relative
- to the bare nucleus in vacuum, ppm
- props.efg -EFG tensors, a.u.^-3
- workaround was implemented, but keep an eye on it.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `textscan()`, `fclose()`, `castep_log()`, `deblank()`, `strcmp()`, `str2double()`, `ischar()`.
