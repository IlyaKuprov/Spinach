# etc/molecules/fatty_acid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/fatty_acid.m`
- Signature: `[sys,inter]=fatty_acid(nprotons)`
- Total lines: 71

## Purpose

Spin system approximating that of a fatty acid. Syntax: [sys,inter]=fatty_acid(nprotons)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(nprotons)`.
- Lines 24-25: Isotope list; implemented by `sys.isotopes=cell(1,nprotons)`.
- Lines 28-29: Chemical shifts; implemented by `inter.zeeman.scalar=cell(1,nprotons)`.
- Lines 34-35: J-couplings; implemented by `inter.coupling.scalar=cell(nprotons,nprotons)`.

### Control flow inferred from the code

- Line 45: `for` loop over `n=4:2:(nprotons-3)`.
- Line 50: `for` loop over `n=5:2:(nprotons-2)`.

### Key state/data transformations

- Lines 25: computes `sys.isotopes` using `sys.isotopes=cell(1,nprotons)`.
- Lines 26: computes `sys.isotopes(:)` using `sys.isotopes(:)={'1H'}`.
- Lines 29: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(1,nprotons)`.
- Lines 30: computes `inter.zeeman.scalar(1:3)` using `inter.zeeman.scalar(1:3)={1.0}`.
- Lines 31: computes `inter.zeeman.scalar(4:2:end)` using `inter.zeeman.scalar(4:2:end)=num2cell(linspace(1.5,5,(nprotons-3)/2))`.
- Lines 32: computes `inter.zeeman.scalar(5:2:end)` using `inter.zeeman.scalar(5:2:end)=num2cell(linspace(1.5,5,(nprotons-3)/2))`.
- Lines 35: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(nprotons,nprotons)`.
- Lines 36: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=20`.
- Lines 37: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=20`.
- Lines 38: computes `inter.coupling.scalar{3,1}` using `inter.coupling.scalar{3,1}=20`.
- Lines 39: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7`.
- Lines 40: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=7`.
- Lines 41: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=7`.
- Lines 42: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=7`.
- Lines 43: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7`.
- Lines 44: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=7`.
- Lines 46: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=15+randn()`.
- Lines 47: computes `inter.coupling.scalar{n,n+2}` using `inter.coupling.scalar{n,n+2}=7 +randn()`.

### Local helper functions

- Line 58: `grumble()` — `function grumble(nprotons)`. Sometimes you see beautiful people with no brains. Sometimes you have ugly people who are intelligent, like scientists.
  - Representative operation: `if (~isnumeric(nprotons))||(~isreal(nprotons))|| (~isscalar(nprotons))||(~isfinite(nprotons))|| (mod(nprotons,1)~=0)||(mod((nprotons-3)/2,1)~=0)|| ((nprotons-3)/2<1)`.
  - Representative operation: `(~isscalar(nprotons))||(~isfinite(nprotons))|| (mod(nprotons,1)~=0)||(mod((nprotons-3)/2,1)~=0)|| ((nprotons-3)/2<1)`.

## Parameters / inputs

- nprotons -the number of protons that the spin
- system should have

## Outputs

- sys, inter -input data structures for Spinach

## Implementation structure

- Spin system approximating that of a fatty acid. Syntax:
- [sys,inter]=fatty_acid(nprotons)
- nprotons -the number of protons that the spin
- system should have
- sys, inter -input data structures for Spinach
- Check consistency
- Isotope list
- Chemical shifts
- J-couplings
- Consistency enforcement
- Sometimes you see beautiful people with no brains. Sometimes
- you have ugly people who are intelligent, like scientists.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `num2cell()`, `isscalar()`.
