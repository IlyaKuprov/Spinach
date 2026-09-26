# kernel/plotting/mri_2d_plot.m

- Signature: `mri_2d_plot(mri,parameters,method)`

## Purpose

2D MRI image plotting. Syntax: mri_2d_plot(mri,parameters,method)

## Physical / mathematical content

## Numerical / algorithmic content

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
