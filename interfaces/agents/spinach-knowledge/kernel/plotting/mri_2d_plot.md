# kernel/plotting/mri_2d_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/mri_2d_plot.m`
- Signature: `mri_2d_plot(mri,parameters,method)`
- Total lines: 185

## Purpose

2D MRI image plotting. Syntax: mri_2d_plot(mri,parameters,method)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(mri,parameters,method)`.
- Lines 41-42: Set the option; implemented by `switch method`.
- Lines 46-47: Get the magnetogyric ratio; implemented by `gamma=spin(parameters.spins{1})`.
- Lines 49-50: Get field of view along dimension 1 (phase encoding); implemented by `max_freq=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp`.
- Lines 53-54: Get field of view along dimension 2 (frequency encoding); implemented by `max_freq=gamma*parameters.ro_grad_dur*parameters.ro_grad_amp`.
- Lines 57-58: Compute axis extents; implemented by `fov1_axis=[-fov1/2,+fov1/2]; fov2_axis=[-fov2/2,+fov2/2]`.
- Lines 60-61: Do the plotting; implemented by `imagesc(fov2_axis,fov1_axis,mri)`.
- Lines 65-66: Get axis extents; implemented by `fov1_axis=[-parameters.dims(1)/2,+parameters.dims(1)/2]`.
- Lines 77-78: Get max spatial frequency along dimension 1 (phase encoding); implemented by `max_freq_1=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp`.
- Lines 80-81: Get max spatial frequency along dimension 2 (frequency encoding); implemented by `max_freq_2=gamma*parameters.ro_grad_dur*parameters.ro_grad_amp`.
- Lines 83-84: Compute axis extents in Hz/m; implemented by `k1_axis=[-max_freq_1/2,+max_freq_1/2]/(2*pi)`.
- Lines 87-88: Plot the real part of the k-space signal; implemented by `imagesc(k2_axis,k1_axis,real(mri))`.
- Lines 92-93: Complain and bomb out; implemented by `error('unknown plotting method.')`.
- Lines 97-98: Do axis labels; implemented by `switch method`.
- Lines 102-103: Fields of view in metres; implemented by `kxlabel('FOV2, m'); kylabel('FOV1, m')`.
- Lines 107-108: Spatial frequency ranges, Hz/m; implemented by `kxlabel('$k_2$, rad/m'); kylabel('$k_1$, rad/m')`.
- Lines 112-113: Apply the colour map; implemented by `colormap(contrast([0 1]))`.
- Lines 115-116: Tidy up axis limits; implemented by `axis equal; xlim tight`.

### Control flow inferred from the code

- Line 42: dispatches on `method`; cases `'image'`, `'phantom'`, `'k-space'`.
- Line 98: dispatches on `method`; cases `{'image','phantom'}`, `'k-space'`.

### Key state/data transformations

- Lines 47: computes `gamma` using `gamma=spin(parameters.spins{1})`.
- Lines 50: computes `max_freq` using `max_freq=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp`.
- Lines 51: computes `pixel_size` using `pixel_size=1*pi/max_freq; fov1=size(mri,1)*pixel_size`.
- Lines 58: computes `fov1_axis` using `fov1_axis=[-fov1/2,+fov1/2]; fov2_axis=[-fov2/2,+fov2/2]`.
- Lines 67: computes `fov2_axis` using `fov2_axis=[-parameters.dims(2)/2,+parameters.dims(2)/2]`.
- Lines 78: computes `max_freq_1` using `max_freq_1=gamma*parameters.pe_grad_dur*parameters.pe_grad_amp`.
- Lines 81: computes `max_freq_2` using `max_freq_2=gamma*parameters.ro_grad_dur*parameters.ro_grad_amp`.
- Lines 84: computes `k1_axis` using `k1_axis=[-max_freq_1/2,+max_freq_1/2]/(2*pi)`.
- Lines 85: computes `k2_axis` using `k2_axis=[-max_freq_2/2,+max_freq_2/2]/(2*pi)`.

### Local helper functions

- Line 122: `grumble()` — `function grumble(mri,parameters,method)`.
  - Representative operation: `if ~ischar(method), error('method must be a character string.'); end`.
  - Representative operation: `if ~ismember(method,{'image','phantom','k-space'})`.

## Parameters / inputs

- method -'image' uses gradient information to
- determine field of view; 'phantom' uses
- phantom dimensions; 'k-space' assumes
- that a k-space representation has been
- supplied and uses gradient information
- and the real part is plotted
- mri -2D MRI image or a phantom, or the k-space
- representation of a 2D MRI image
- parameters.spins -nuclei on which the sequence ran.
- parameters.pe_grad_amp -amplitude of the phase encoding
- gradient, T/m
- parameters.pe_grad_dur -duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_amp -the amplitude of the readout
- gradient, T/m
- parameters.ro_grad_dur -the duration of the readout
- gradient, seconds.

## Implementation structure

- 2D MRI image plotting. Syntax:
- mri_2d_plot(mri,parameters,method)
- method -'image' uses gradient information to
- determine field of view; 'phantom' uses
- phantom dimensions; 'k-space' assumes
- that a k-space representation has been
- supplied and uses gradient information
- and the real part is plotted
- mri -2D MRI image or a phantom, or the k-space
- representation of a 2D MRI image
- parameters.spins -nuclei on which the sequence ran.
- parameters.pe_grad_amp -amplitude of the phase encoding

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `kxlabel()`, `kylabel()`, `colormap()`, `contrast()`, `ischar()`, `ismember()`, `ismatrix()`, `strcmp()`, `isfield()`, `iscell()`.
