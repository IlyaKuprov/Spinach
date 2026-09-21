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
