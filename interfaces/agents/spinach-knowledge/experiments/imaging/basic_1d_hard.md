# experiments/imaging/basic_1d_hard.m

- Signature: `fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)`

## Purpose

Basic 1D imaging sequence with a hard pulse. Syntax: fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.ro_grad_amp -readout gradient amplitude, T/m parameters.sweep -detection sweep width, Hz parameters.npoints -number of points in the fid parameters.offset -transmitter and receiver offset, Hz

## Physical / mathematical content

A hard 90-degree pulse about y is followed by evolution, a hard 180-degree pulse about x, pre-phasing with `-ro_grad_amp*G{1}`, and acquisition with `+ro_grad_amp*G{1}`.

## Numerical / algorithmic content

The function forms `L=H+F+1i*R+1i*K`, applies the specified pulses and evolutions, and returns the acquired FID; it does not perform an FFT.

## Outputs

- fid -free induction decay that should be Fourier transformed
- to obtain the image

## Implementation structure

The local input checks validate the formalism, operator dimensions, gradient cell, and required parameters. The function builds the pulse operators, applies the sequence, and returns `acquire` output.
