# etc/molecules/allyl_pyruvate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/allyl_pyruvate.m`
- Signature: `[sys,inter]=allyl_pyruvate(spins)`
- Total lines: 165

## Purpose

Spin system of allyl pyruvate. Isotropic chemical shifts and J-couplings determined by spectral fitting, coordinates and chemical shift anisotropies from DFT. Syntax: [sys,inter]=allyl_pyruvate(spins)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(spins)`.
- Lines 30-32: Spin system; implemented by `sys.isotopes={'13C','13C','13C','13C','13C','13C', '1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 36-39: Chemical shifts (from fitting); implemented by `inter.zeeman.scalar={119.9643 130.6962 66.8463 160.4202 191.6305 26.6961 5.4086 5.3273 5.9678 4.7466 4.7466 2.4836 2.4836 2.4836}`.
- Lines 41-42: Chemical shift anisotropies (DFT); implemented by `inter.zeeman.matrix{idxof(sys,'C1')}=-remtrace([ 73.8088 56.3125 -5.7848`.
- Lines 85-86: 13C-1H J-couplings (from fitting); implemented by `inter.coupling.scalar=cell(14,14)`.
- Lines 113-114: 1H-1H J-couplings (from fitting); implemented by `inter.coupling.scalar{idxof(sys,'Ha'),idxof(sys,'Hb')}= +1.16`.
- Lines 124-125: Coordinates (from DFT); implemented by `inter.coordinates={[-3.845334 -0.310277 -0.426714]`.
- Lines 140-141: Decide which atoms to return; implemented by `mask=ismember(sys.isotopes,spins)`.
- Lines 143-144: Prune the arrays; implemented by `sys.isotopes=sys.isotopes(mask)`.

### Key state/data transformations

- Lines 31-32: computes `sys.isotopes` using `sys.isotopes={'13C','13C','13C','13C','13C','13C', '1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 33-34: computes `sys.labels` using `sys.labels={'C1','C2','C3','C4','C5','C6', 'Ha','Hb','Hc','Hd1','Hd2','He1','He2','He3'}`.
- Lines 37-39: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={119.9643 130.6962 66.8463 160.4202 191.6305 26.6961 5.4086 5.3273 5.9678 4.7466 4.7466 2.4836 2.4836 2.4836}`.
- Lines 42: computes `inter.zeeman.matrix{idxof(sys,'C1')}` using `inter.zeeman.matrix{idxof(sys,'C1')}=-remtrace([ 73.8088 56.3125 -5.7848`.
- Lines 45: computes `inter.zeeman.matrix{idxof(sys,'C2')}` using `inter.zeeman.matrix{idxof(sys,'C2')}=-remtrace([ 65.5131 69.8026 0.8618`.
- Lines 48: computes `inter.zeeman.matrix{idxof(sys,'C3')}` using `inter.zeeman.matrix{idxof(sys,'C3')}=-remtrace([140.4211 17.8082 1.7258`.
- Lines 51: computes `inter.zeeman.matrix{idxof(sys,'C4')}` using `inter.zeeman.matrix{idxof(sys,'C4')}=-remtrace([-64.9706 -30.5508 -9.0450`.
- Lines 54: computes `inter.zeeman.matrix{idxof(sys,'C5')}` using `inter.zeeman.matrix{idxof(sys,'C5')}=-remtrace([-83.4746 -39.6016 -10.1866`.
- Lines 57: computes `inter.zeeman.matrix{idxof(sys,'C6')}` using `inter.zeeman.matrix{idxof(sys,'C6')}=-remtrace([184.8168 -3.1537 1.8807`.
- Lines 60: computes `inter.zeeman.matrix{idxof(sys,'Ha')}` using `inter.zeeman.matrix{idxof(sys,'Ha')}=-remtrace([ 28.2575 1.3410 -3.0000`.
- Lines 63: computes `inter.zeeman.matrix{idxof(sys,'Hb')}` using `inter.zeeman.matrix{idxof(sys,'Hb')}=-remtrace([ 28.0342 -0.1603 -1.2955`.
- Lines 66: computes `inter.zeeman.matrix{idxof(sys,'Hc')}` using `inter.zeeman.matrix{idxof(sys,'Hc')}=-remtrace([ 27.9364 -0.1637 -1.6581`.
- Lines 69: computes `inter.zeeman.matrix{idxof(sys,'Hd1')}` using `inter.zeeman.matrix{idxof(sys,'Hd1')}=-remtrace([28.0196 1.6136 1.4936`.
- Lines 72: computes `inter.zeeman.matrix{idxof(sys,'Hd2')}` using `inter.zeeman.matrix{idxof(sys,'Hd2')}=-remtrace([29.6956 1.5093 -1.5351`.
- Lines 75: computes `inter.zeeman.matrix{idxof(sys,'He1')}` using `inter.zeeman.matrix{idxof(sys,'He1')}=-remtrace([29.3456 -0.8213 -1.3863`.
- Lines 78: computes `inter.zeeman.matrix{idxof(sys,'He2')}` using `inter.zeeman.matrix{idxof(sys,'He2')}=-remtrace([28.9925 -0.4100 1.3818`.
- Lines 81: computes `inter.zeeman.matrix{idxof(sys,'He3')}` using `inter.zeeman.matrix{idxof(sys,'He3')}=-remtrace([32.9096 4.0247 0.3478`.
- Lines 86: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(14,14)`.

### Local helper functions

- Line 154: `grumble()` — `function grumble(spins)`. You will sag in your chains as embers scorch your feet. It will be endless. The word 'endless' -do
  - Representative operation: `if (~iscell(spins))||(~all(cellfun(@ischar,spins)))`.
  - Representative operation: `error('spins argument must be a cell array of character strings.')`.

## Parameters / inputs

- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}

## Outputs

- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.

## Implementation structure

- Spin system of allyl pyruvate. Isotropic chemical shifts and
- J-couplings determined by spectral fitting, coordinates and
- chemical shift anisotropies from DFT. Syntax:
- [sys,inter]=allyl_pyruvate(spins)
- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}
- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.
- Check consistency
- Spin system

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `idxof()`, `remtrace()`, `ismember()`, `iscell()`, `all()`, `cellfun()`.
