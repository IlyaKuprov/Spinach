# kernel/utilities/idxof.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/idxof.m`
- Signature: `idx=idxof(sys,label)`
- Total lines: 49

## Purpose

Allows interaction specification by spin label rather than number. Syntax: idx=idxof(sys,label)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(sys,label)`.
- Lines 27-28: Locate the label; implemented by `idx=find(cellfun(@(x)strcmp(label,x),sys.labels))`.
- Lines 30-31: Check the output; implemented by `if isempty(idx), error('label not found.'); end`.

### Control flow inferred from the code

- Line 31: conditional branch on `isempty(idx), error('label not found.'); end`.
- Line 32: conditional branch on `numel(idx)>1, error('labels are not unique'); end`.

### Key state/data transformations

- Lines 28: computes `idx` using `idx=find(cellfun(@(x)strcmp(label,x),sys.labels))`.

### Local helper functions

- Line 37: `grumble()` — `function grumble(sys,label)`. Советскую власть я возненавидела ещё до того, как узнала, что она есть.
  - Representative operation: `if ~ischar(label), error('label must be a character string.'); end`.
  - Representative operation: `if (~isfield(sys,'labels')), error('sys.labels is missing.'); end`.

## Parameters / inputs

- sys -Spinach input structure that
- includes a sys.labels field
- with unique labels
- label -label whose index is to be returned

## Outputs

- idx -the index of the spin, an integer

## Implementation structure

- Allows interaction specification by spin label
- rather than number. Syntax:
- idx=idxof(sys,label)
- sys -Spinach input structure that
- includes a sys.labels field
- with unique labels
- label -label whose index is to be returned
- idx -the index of the spin, an integer
- Check consistency
- Locate the label
- Check the output
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `strcmp()`, `ischar()`, `isfield()`.
