# interfaces/spinxml/x2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/spinxml/x2spinach.m`
- Signature: `[sys,inter]=x2spinach(filename,shielding_refs)`
- Total lines: 895

## Purpose

Reads SpinXML files and forms Spinach data structures. Syntax: [sys,inter]=x2spinach(file_name,shielding_refs)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(filename,shielding_refs)`.
- Lines 33-34: Get the data structures started; implemented by `spin_ids=[]; sys.isotopes={}`.
- Lines 37-38: Parse the XML; implemented by `xml=parsexml(filename)`.
- Lines 40-41: First pass -discover spins; implemented by `for n=1:numel(xml.children)`.
- Lines 43-44: Look for spin elements; implemented by `if strcmpi(xml.children(n).name,'spin')`.
- Lines 46-47: Labels and coordinates are optional; implemented by `current_spin_label=''; xyz=[]`.
- Lines 49-50: Look for spin attributes; implemented by `for k=1:numel(xml.children(n).attributes)`.
- Lines 52-53: Process attributes; implemented by `if strcmpi(xml.children(n).attributes(k).name,'id')`.
- Lines 55-56: Process spin id number; implemented by `current_spin_id=str2double(xml.children(n).attributes(k).value)`.
- Lines 58-59: Check for zero-base indexing; implemented by `if current_spin_id==0`.
- Lines 65-66: Process isotope specifications; implemented by `current_spin_isotope=xml.children(n).attributes(k).value`.
- Lines 70-71: Process text labels; implemented by `current_spin_label=xml.children(n).attributes(k).value`.
- Lines 77-78: Look for coordinate information; implemented by `for k=1:numel(xml.children(n).children)`.
- Lines 80-81: Parse coordinate information; implemented by `if strcmpi(xml.children(n).children(k).name,'coordinates')`.
- Lines 83-84: Process Cartesian vector components; implemented by `for m=1:numel(xml.children(n).children(k).attributes)`.
- Lines 98-99: Form Spinach structures; implemented by `spin_ids=[spin_ids current_spin_id]`.
- Lines 104-105: Clear temporary variables; implemented by `clear('xyz','current_spin_id','current_spin_isotope','current_spin_label')`.
- Lines 111-112: Sort arrays by spin index number; implemented by `[spin_ids,index]=sort(spin_ids,'ascend')`.

### Control flow inferred from the code

- Line 41: `for` loop over `n=1:numel(xml.children)`.
- Line 44: conditional branch on `strcmpi(xml.children(n).name,'spin')`.
- Line 50: `for` loop over `k=1:numel(xml.children(n).attributes)`.
- Line 53: conditional branch on `strcmpi(xml.children(n).attributes(k).name,'id')`.
- Line 59: conditional branch on `current_spin_id==0`.
- Line 78: `for` loop over `k=1:numel(xml.children(n).children)`.
- Line 81: conditional branch on `strcmpi(xml.children(n).children(k).name,'coordinates')`.
- Line 84: `for` loop over `m=1:numel(xml.children(n).children(k).attributes)`.
- Line 85: conditional branch on `strcmpi(xml.children(n).children(k).attributes(m).name,'x')`.
- Line 118: conditional branch on `~isequal(spin_ids(:).',1:numel(spin_ids))`.
- Line 135: `for` loop over `n=1:numel(xml.children)`.
- Line 138: conditional branch on `strcmpi(xml.children(n).name,'interaction')`.
- Line 145: `for` loop over `k=1:numel(xml.children(n).attributes)`.
- Line 148: conditional branch on `strcmpi(xml.children(n).attributes(k).name,'kind')`.

### Key state/data transformations

- Lines 34: computes `spin_ids` using `spin_ids=[]; sys.isotopes={}`.
- Lines 35: computes `sys.labels` using `sys.labels={}; inter.coordinates={}`.
- Lines 38: computes `xml` using `xml=parsexml(filename)`.
- Lines 47: computes `current_spin_label` using `current_spin_label=''; xyz=[]`.
- Lines 56: computes `current_spin_id` using `current_spin_id=str2double(xml.children(n).attributes(k).value)`.
- Lines 66: computes `current_spin_isotope` using `current_spin_isotope=xml.children(n).attributes(k).value`.
- Lines 86: computes `xyz(1)` using `xyz(1)=str2double(xml.children(n).children(k).attributes(m).value)`.
- Lines 88: computes `xyz(2)` using `xyz(2)=str2double(xml.children(n).children(k).attributes(m).value)`.
- Lines 90: computes `xyz(3)` using `xyz(3)=str2double(xml.children(n).children(k).attributes(m).value)`.
- Lines 100: computes `sys.isotopes` using `sys.isotopes=[sys.isotopes {current_spin_isotope}]`.
- Lines 102: computes `inter.coordinates` using `inter.coordinates=[inter.coordinates; {xyz}]`.
- Lines 112: computes `[spin_ids,index]` using `[spin_ids,index]=sort(spin_ids,'ascend')`.
- Lines 123: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 124: computes `inter.zeeman.reference` using `inter.zeeman.reference=cell(1,numel(sys.isotopes))`.
- Lines 125: computes `inter.zeeman.label` using `inter.zeeman.label=cell(1,numel(sys.isotopes))`.
- Lines 126: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 127: computes `inter.coupling.label` using `inter.coupling.label=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 128: computes `inter.spinrot.matrix` using `inter.spinrot.matrix=cell(1,numel(sys.isotopes))`.

### Local helper functions

- Line 865: `grumble()` — `function grumble(filename,shielding_refs)`.
  - Representative operation: `if (~ischar(filename))||isempty(filename)`.
  - Representative operation: `error('filename must be a non-empty character string.')`.

## Parameters / inputs

- shielding_refs -when absolute shielding is provided, this
- array must give absolute shielding of the
- corresponding nuclei in the reference sub-
- stance, e.g. {{'1H',31.5},{'13C',189.7}}.
- This is necessary because Spinach requires
- chemical shifts rather than shieldings.
- If chemical shifts are provided, use an
- empty cell array.

## Outputs

- sys, inter -Spinach spin system input structures
- WARNING: this function assumes that the SpinXML file has passed
- the validation against the schema, which may be obtain-
- ed from http://spindynamics.org/SpinXML.php

## Implementation structure

- Reads SpinXML files and forms Spinach data structures. Syntax:
- [sys,inter]=x2spinach(file_name,shielding_refs)
- shielding_refs -when absolute shielding is provided, this
- array must give absolute shielding of the
- corresponding nuclei in the reference sub-
- stance, e.g. {{'1H',31.5},{'13C',189.7}}.
- This is necessary because Spinach requires
- chemical shifts rather than shieldings.
- If chemical shifts are provided, use an
- empty cell array.
- sys, inter -Spinach spin system input structures
- WARNING: this function assumes that the SpinXML file has passed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `parsexml()`, `strcmpi()`, `str2double()`, `xyz()`, `clear()`, `isequal()`, `spin_ids()`, `anas2mat()`, `axrh2mat()`, `spsk2mat()`, `euler2dcm()`, `qter2dcm()`, `dcm()`, `anax2dcm()`, `iselectron()`.
