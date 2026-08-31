# interfaces/gissmo/gissmo2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gissmo/gissmo2spinach.m`
- Signature: `[sys,inter]=gissmo2spinach(filename,subsystem)`
- Total lines: 206

## Purpose

Reads GISSMO files and forms Spinach data structures. Syntax: [sys,inter]=gissmo2spinach(file_name,subsystem)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(filename)`.
- Lines 32-33: Parse the XML; implemented by `disp('parsing GISSMO xml ')`.
- Lines 36-37: Set essential flags; implemented by `magnet=false(); lw=false()`.
- Lines 41-42: Look inside; implemented by `for n=1:numel(xml.children)`.
- Lines 44-45: Look for magnetic induction; implemented by `if strcmpi(xml.children(n).name,'field_strength')`.
- Lines 47-48: Get magnetic induction; implemented by `sys.magnet=2*pi*1e6*str2double(xml.children(n).children.data)/spin('1H')`.
- Lines 50-51: Inform the user; implemented by `disp('found magnet induction'); magnet=true()`.
- Lines 55-57: Look for a coupling matrix; implemented by `if strcmpi(xml.children(n).name,'coupling_matrix')&& (current_subsystem==subsystem)`.
- Lines 59-60: Grab a shorthand; implemented by `cm=xml.children(n)`.
- Lines 62-63: Look inside; implemented by `for k=1:numel(cm.children)`.
- Lines 65-66: Look for line width; implemented by `if strcmpi(cm.children(k).name,'lw')`.
- Lines 68-69: Get the line width; implemented by `inter.relaxation={'damp'}`.
- Lines 74-75: Inform the user; implemented by `disp('found line width'); lw=true()`.
- Lines 79-80: Look for the spin list; implemented by `if strcmpi(cm.children(k).name,'spin_names')`.
- Lines 82-83: Grab a shorthand; implemented by `sn=cm.children(k)`.
- Lines 85-86: Look inside; implemented by `for m=1:numel(sn.children)`.
- Lines 88-89: Look for spin names; implemented by `if strcmpi(sn.children(m).name,'spin')`.
- Lines 91-92: Assign labels; implemented by `idx=strcmpi('index',{sn.children(m).attributes.name})`.

### Control flow inferred from the code

- Line 42: `for` loop over `n=1:numel(xml.children)`.
- Line 45: conditional branch on `strcmpi(xml.children(n).name,'field_strength')`.
- Line 56: conditional branch on `strcmpi(xml.children(n).name,'coupling_matrix')&&`.
- Line 63: `for` loop over `k=1:numel(cm.children)`.
- Line 66: conditional branch on `strcmpi(cm.children(k).name,'lw')`.
- Line 80: conditional branch on `strcmpi(cm.children(k).name,'spin_names')`.
- Line 86: `for` loop over `m=1:numel(sn.children)`.
- Line 89: conditional branch on `strcmpi(sn.children(m).name,'spin')`.
- Line 108: conditional branch on `strcmpi(cm.children(k).name,'chemical_shifts_ppm')`.
- Line 114: `for` loop over `m=1:numel(cs.children)`.
- Line 117: conditional branch on `strcmpi(cs.children(m).name,'cs')`.
- Line 136: conditional branch on `strcmpi(cm.children(k).name,'couplings_hz')`.
- Line 146: `for` loop over `m=1:numel(jc.children)`.
- Line 149: conditional branch on `strcmpi(jc.children(m).name,'coupling')`.

### Key state/data transformations

- Lines 34: computes `xml` using `xml=parsexml(filename)`.
- Lines 37: computes `magnet` using `magnet=false(); lw=false()`.
- Lines 38: computes `shifts` using `shifts=false(); jc=false()`.
- Lines 39: computes `current_subsystem` using `current_subsystem=1`.
- Lines 48: computes `sys.magnet` using `sys.magnet=2*pi*1e6*str2double(xml.children(n).children.data)/spin('1H')`.
- Lines 51: computes `disp('found magnet induction'); magnet` using `disp('found magnet induction'); magnet=true()`.
- Lines 60: computes `cm` using `cm=xml.children(n)`.
- Lines 69: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 70: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 71: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 72: computes `inter.damp_rate` using `inter.damp_rate=fwhm2rlx(str2double(cm.children(k).children.data))`.
- Lines 75: computes `disp('found line width'); lw` using `disp('found line width'); lw=true()`.
- Lines 83: computes `sn` using `sn=cm.children(k)`.
- Lines 92: computes `idx` using `idx=strcmpi('index',{sn.children(m).attributes.name})`.
- Lines 93: computes `lab` using `lab=strcmpi('name',{sn.children(m).attributes.name})`.
- Lines 96: computes `sys.labels{idx}` using `sys.labels{idx}=['Atom ' lab]`.
- Lines 111: computes `cs` using `cs=cm.children(k)`.
- Lines 121: computes `ppm` using `ppm=strcmpi('ppm',{cs.children(m).attributes.name})`.

### Local helper functions

- Line 192: `grumble()` — `function grumble(filename)`. Perhaps this habit goes back to the primitive belief that the word and the thing, the name and the object,
  - Representative operation: `if (~ischar(filename))||isempty(filename)`.
  - Representative operation: `error('filename must be a non-empty character string.')`.

## Parameters / inputs

- file_name -character string with the name of
- the GISSMO XML file
- subsystem -which of the coupling matrices to
- to import

## Outputs

- sys, inter -Spinach data structures, ready for
- calling create.m
- Note: GISSMO only provides chemical shifts, J-couplings, the
- non-selective line width, and the magnet field. You may
- want to add further parameters by editing sys and inter
- data structures manually.

## Implementation structure

- Reads GISSMO files and forms Spinach data structures. Syntax:
- [sys,inter]=gissmo2spinach(file_name,subsystem)
- file_name - character string with the name of
- the GISSMO XML file
- subsystem - which of the coupling matrices to
- to import
- sys, inter - Spinach data structures, ready for
- calling create.m
- Note: GISSMO only provides chemical shifts, J-couplings, the
- non-selective line width, and the magnet field. You may
- want to add further parameters by editing sys and inter
- data structures manually.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `parsexml()`, `false()`, `strcmpi()`, `str2double()`, `spin()`, `true()`, `fwhm2rlx()`, `ischar()`, `exist()`.
