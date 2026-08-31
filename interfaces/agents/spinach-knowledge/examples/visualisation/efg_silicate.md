# examples/visualisation/efg_silicate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/visualisation/efg_silicate.m`
- Signature: `efg_silicate()`
- Total lines: 25

## Purpose

Example of electric field gradient tensor visualisation for an aluminosilicate solid. CASTEP log is parsed.

## Physical / mathematical content

- Visualisation examples. These scripts expose tensor geometry: principal axes, shielding/hyperfine/EFG ellipsoids, molecular frames, and the relationship between tensor eigenstructure and observable anisotropy.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Import CASTEP data; implemented by `props=c2spinach('alsilicate.magres')`.
- Lines 13-14: Do the visualisation; implemented by `kfigure(); scale_figure([1.875 1.125])`.

### Key state/data transformations

- Lines 11: computes `props` using `props=c2spinach('alsilicate.magres')`.
- Lines 15: computes `subplot(1,2,1); options.style` using `subplot(1,2,1); options.style='ellipsoids'`.
- Lines 19: computes `subplot(1,2,2); options.style` using `subplot(1,2,2); options.style='harmonics'`.

## Implementation structure

- Example of electric field gradient tensor visualisation for
- an aluminosilicate solid. CASTEP log is parsed.
- Import CASTEP data
- Do the visualisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `c2spinach()`, `kfigure()`, `scale_figure()`, `subplot()`, `efg_display()`, `set()`, `ktitle()`.
