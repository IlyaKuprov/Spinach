# kernel/plotting/axis_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/axis_1d.m`
- Signature: `[ax,ax_label]=axis_1d(spin_system,parameters)`
- Total lines: 220

## Purpose

Generates axis ticks for plotting 1D spectra. Syntax: [ax,ax_label]=axis_1d(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 41-42: Build the axis; implemented by `if isscalar(parameters.sweep)`.
- Lines 44-45: Offset + sweep specification; implemented by `ax=ft_axis(parameters.offset,parameters.sweep,parameters.zerofill)`.
- Lines 47-48: Report to the user; implemented by `report(spin_system,'offset + sweep frequency axis specification:')`.
- Lines 54-55: Two-element sweep specification; implemented by `ax=linspace(parameters.sweep(1),parameters.sweep(2),parameters.zerofill)`.
- Lines 57-58: Report to the user; implemented by `report(spin_system,'two-element sweep frequency axis specification:')`.
- Lines 64-65: Carrier frequency of the spin in Hz; implemented by `basefrq=-spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 67-68: Make the isotope mass number an upper index; implemented by `k=regexp(parameters.spins{1},'[A-Za-z]','once')`.
- Lines 76-77: Convert the units if necessary; implemented by `switch parameters.axis_units`.
- Lines 81-82: NMR-style ppm scale with respect to the carrier frequency; implemented by `ax_label=[iso_label ' chemical shift / ppm']`.
- Lines 87-88: Pulsed FFT EPR style Gauss axis for rotating frame simulations; implemented by `ax_label='Magnetic induction / Gauss'`.
- Lines 93-94: Pulsed FFT EPR style mT axis for rotating frame simulations; implemented by `ax_label='Magnetic induction / mT'`.
- Lines 99-100: Raw Hz axis for rotating frame simulations; implemented by `ax_label=[iso_label ' offset frequency / Hz']`.
- Lines 104-105: Raw kHz axis for rotating frame simulations; implemented by `ax_label=[iso_label ' offset frequency / kHz']`.
- Lines 110-111: Raw MHz axis for rotating frame simulations; implemented by `ax_label=[iso_label ' offset frequency / MHz']`.
- Lines 116-117: Raw MHz axis for lab frame simulations; implemented by `ax_label=[iso_label ' frequency / MHz']`.
- Lines 122-123: Raw GHz axis for rotating frame simulations; implemented by `ax_label=[iso_label ' offset frequency / GHz']`.
- Lines 128-129: Raw GHz axis for lab frame simulations; implemented by `ax_label=[iso_label ' frequency / GHz']`.

### Control flow inferred from the code

- Line 42: conditional branch on `isscalar(parameters.sweep)`.
- Line 69: conditional branch on `k>1`.
- Line 77: dispatches on `parameters.axis_units`; cases `'ppm'`, `'Gauss'`, `'mT'`, `'Hz'`, `'kHz'`, `'MHz'`, `'MHz-labframe'`, `'GHz'`, `'GHz-labframe'`, `'gtensor'`, ….

### Key state/data transformations

- Lines 45: computes `ax` using `ax=ft_axis(parameters.offset,parameters.sweep,parameters.zerofill)`.
- Lines 65: computes `basefrq` using `basefrq=-spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`.
- Lines 68: computes `k` using `k=regexp(parameters.spins{1},'[A-Za-z]','once')`.
- Lines 70-71: computes `iso_label` using `iso_label=sprintf('$^{%s}$%s',parameters.spins{1}(1:k-1), parameters.spins{1}(k:end))`.
- Lines 82: computes `ax_label` using `ax_label=[iso_label ' chemical shift / ppm']`.

### Local helper functions

- Line 154: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'sweep')`.
  - Representative operation: `error('sweep width should be specified in parameters.sweep variable.')`.

## Parameters / inputs

- parameters.sweep -either a one-element array giving the sweep width
- in Hz, or a two-element array giving the spectral
- extents in Hz around the offset.
- parameters.zerofill -the number of points in the NMR spectrum after
- zerofilling and Fourier transform
- parameters.offset -offset of the spectrum centre point relative to
- the magnet frequency, Hz
- parameters.axis_units -a character string with the units in which the
- axis ticks should be returned: 'ppm', 'Gauss',
- 'mT', 'Hz', 'kHz', 'MHz','MHz-labframe', 'GHz',
- 'GHz-labframe', 'gtensor', 'points'
- parameters.spins -the spin involved, e.g. {'1H'}

## Outputs

- axis -a row vector of axis tick values
- ax_label -axis label for displaying on the plot
- Note: magnetic field units use the free electron g-tensor for conversion.

## Implementation structure

- Generates axis ticks for plotting 1D spectra. Syntax:
- [ax,ax_label]=axis_1d(spin_system,parameters)
- parameters.sweep - either a one-element array giving the sweep width
- in Hz, or a two-element array giving the spectral
- extents in Hz around the offset.
- parameters.zerofill -the number of points in the NMR spectrum after
- zerofilling and Fourier transform
- parameters.offset -offset of the spectrum centre point relative to
- the magnet frequency, Hz
- parameters.axis_units -a character string with the units in which the
- axis ticks should be returned: 'ppm', 'Gauss',
- 'mT', 'Hz', 'kHz', 'MHz','MHz-labframe', 'GHz',

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `ft_axis()`, `report()`, `num2str()`, `spin()`, `regexp()`, `isfield()`, `ischar()`, `iscell()`, `ismember()`, `strcmp()`.
