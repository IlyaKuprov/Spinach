# kernel/pulses/read_wave.m

- Signature: `[A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)`

## Purpose

Reads JCAMP-DX pulse waveform files (a few examples are distri- buted with Spinach, see /kernel/pulses/pk_files). Syntax: [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- filename -a string containing the name of the file
- npoints -waveform upsampling or downsampling is
- performed to this number of points

## Outputs

- A -polar amplitude at each slice
- phi -polar phase at each slice, radians
- Cx -Cartesian amplitude in X at each slice
- Cy -Cartesian amplitude in Y at each slice
- scaling factor -scaling factor for a given pulse shape
- Note: put your own pulses into /kernel/pulses/pk_files; please
- also consider sending them to us.

## Implementation structure

- Reads JCAMP-DX pulse waveform files (a few examples are distri-
- buted with Spinach, see /kernel/pulses/pk_files). Syntax:
- [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)
- filename - a string containing the name of the file
- npoints - waveform upsampling or downsampling is
- performed to this number of points
- A - polar amplitude at each slice
- phi - polar phase at each slice, radians
- Cx - Cartesian amplitude in X at each slice
- Cy - Cartesian amplitude in Y at each slice
- scaling factor - scaling factor for a given pulse shape
- Note: put your own pulses into /kernel/pulses/pk_files; please
