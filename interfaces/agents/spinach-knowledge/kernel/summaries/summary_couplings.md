# kernel/summaries/summary_couplings.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_couplings.m`
- Signature: `summary_couplings(spin_system,header)`
- Total lines: 82

## Purpose

Prints spin-spin coupling tensor summary for a Spinach system. Syntax: summary_couplings(spin_system,header)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(spin_system,header)`.
- Lines 25-26: Print the summary table; implemented by `report(spin_system,header)`.
- Lines 31-32: Detect significant couplings; implemented by `[rows,cols,~]=find(cellfun(@(x)norm(x,2),spin_system.inter.coupling.matrix)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 34-35: Loop over significant couplings; implemented by `for n=1:numel(rows)`.
- Lines 37-38: Get the isotropic part; implemented by `iso=trace(spin_system.inter.coupling.matrix{rows(n),cols(n)})/3`.
- Lines 40-41: Get the first and second rank parts; implemented by `[~,rank1,rank2]=mat2sphten(spin_system.inter.coupling.matrix{rows(n),cols(n)})`.
- Lines 44-47: Do the printing in Hz; implemented by `report(spin_system,[pad(num2str(rows(n)),8) pad(num2str(cols(n)),8) num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(1,:)/(2*pi),'%+10.6e %+10.6e %+10.6e')])`.
- Lines 57-58: Print the break line; implemented by `if n<numel(rows), report(spin_system,' '); end`.

### Control flow inferred from the code

- Line 35: `for` loop over `n=1:numel(rows)`.
- Line 58: conditional branch on `n<numel(rows), report(spin_system,' '); end`.

### Key state/data transformations

- Lines 32: computes `[rows,cols,~]` using `[rows,cols,~]=find(cellfun(@(x)norm(x,2),spin_system.inter.coupling.matrix)>2*pi*spin_system.tols.inter_cutoff)`.
- Lines 38: computes `iso` using `iso=trace(spin_system.inter.coupling.matrix{rows(n),cols(n)})/3`.
- Lines 41: computes `[~,rank1,rank2]` using `[~,rank1,rank2]=mat2sphten(spin_system.inter.coupling.matrix{rows(n),cols(n)})`.
- Lines 42: computes `rank1` using `rank1=sphten2mat([],rank1,[]); rank2=sphten2mat([],[],rank2)`.

### Local helper functions

- Line 66: `grumble()` — `function grumble(spin_system,header)`. A gentleman a few rows in front of us took grave exception to the behaviour of an opposing player
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints spin-spin coupling tensor summary for a Spinach system. Syntax:
- summary_couplings(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Detect significant couplings
- Loop over significant couplings
- Get the isotropic part
- Get the first and second rank parts

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `cellfun()`, `rows()`, `cols()`, `mat2sphten()`, `sphten2mat()`, `pad()`, `num2str()`, `isstruct()`, `ischar()`.
