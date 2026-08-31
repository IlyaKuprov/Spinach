# examples/visualisation/hfc_porphyrine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/visualisation/hfc_porphyrine.m`
- Signature: `hfc_porphyrine()`
- Total lines: 26

## Purpose

Example of proton hyperfine tensor visualisation for copper porphyrine. ORCA log is parsed.

## Physical / mathematical content

- Visualisation examples. These scripts expose tensor geometry: principal axes, shielding/hyperfine/EFG ellipsoids, molecular frames, and the relationship between tensor eigenstructure and observable anisotropy.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Read the ORCA log; implemented by `props=oparse('porphyrine.out')`.
- Lines 12-13: Do the visualisation; implemented by `kfigure(); subplot(1,2,1)`.

### Key state/data transformations

- Lines 10: computes `props` using `props=oparse('porphyrine.out')`.
- Lines 14: computes `options.style` using `options.style='ellipsoids'`.

## Implementation structure

- Example of proton hyperfine tensor visualisation for
- copper porphyrine. ORCA log is parsed.
- Read the ORCA log
- Do the visualisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `oparse()`, `kfigure()`, `subplot()`, `hfc_display()`, `set()`, `ktitle()`, `scale_figure()`.
