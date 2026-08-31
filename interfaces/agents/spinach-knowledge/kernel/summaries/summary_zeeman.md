# kernel/summaries/summary_zeeman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_zeeman.m`
- Signature: `summary_zeeman(spin_system,header)`
- Total lines: 78

## Purpose

Prints Zeeman interaction tensor summary for a Spinach system. Syntax: summary_zeeman(spin_system,header)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(spin_system,header)`.
- Lines 25-26: Print the summary table; implemented by `report(spin_system,header)`.
- Lines 33-34: Get the isotropic part; implemented by `iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 36-37: Get the first and second rank parts; implemented by `[~,rank1,rank2]=mat2sphten(spin_system.inter.zeeman.matrix{n})`.
- Lines 40-44: Do the printing; implemented by `report(spin_system,[pad(num2str(n),5) pad(spin_system.comp.isotopes{n},6) pad(num2str(spin_system.comp.mults(n)),7) num2str(spin_system.inter.zeeman.matrix{n}(1,:),'%+10…`.
- Lines 54-55: Print the break line; implemented by `if n<spin_system.comp.nspins, report(spin_system,' '); end`.

### Control flow inferred from the code

- Line 30: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 31: conditional branch on `~isempty(spin_system.inter.zeeman.matrix{n})`.
- Line 55: conditional branch on `n<spin_system.comp.nspins, report(spin_system,' '); end`.

### Key state/data transformations

- Lines 34: computes `iso` using `iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
- Lines 37: computes `[~,rank1,rank2]` using `[~,rank1,rank2]=mat2sphten(spin_system.inter.zeeman.matrix{n})`.
- Lines 38: computes `rank1` using `rank1=sphten2mat([],rank1,[]); rank2=sphten2mat([],[],rank2)`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(spin_system,header)`. Jim wants me to go out with other men so that he will have something to write about.
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints Zeeman interaction tensor summary for a Spinach system. Syntax:
- summary_zeeman(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Get the isotropic part
- Get the first and second rank parts
- Do the printing
- Print the break line

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `mat2sphten()`, `sphten2mat()`, `pad()`, `num2str()`, `isstruct()`, `ischar()`.
