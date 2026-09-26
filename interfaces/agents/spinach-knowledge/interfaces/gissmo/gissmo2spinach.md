# interfaces/gissmo/gissmo2spinach.m

- Signature: `[sys,inter]=gissmo2spinach(filename,subsystem)`

## Purpose

Reads a GISSMO XML file and forms the Spinach system and interaction structures for the selected coupling matrix.

## Physical / mathematical content

GISSMO provides proton chemical shifts in ppm, scalar J-couplings in hertz, the proton spectrometer frequency, and a non-selective linewidth. The frequency determines the magnetic induction. The linewidth is treated as Lorentzian FWHM in hertz and converted by `fwhm2rlx` to the damping rate `pi*FWHM` in inverse seconds.

Pure non-selective damping uses zero equilibrium and full (`labframe`) retention. In either Liouville formalism its generator is `-rate*(I-u*u')`, where `u` is the normalised unit state: every traceless state decays, while trace and the identity are preserved. Damping is added after retention, so this setting preserves the previous spherical-tensor result without requesting unsupported Zeeman diagonal retention; it does not require a laboratory-frame Hamiltonian for pure damping.

## Numerical / algorithmic content

The XML coupling matrices are selected by one-based order. Spin labels, shifts, and couplings use their XML indices; the imported spins are protons. Essential magnetic-field, linewidth, shift, and coupling data must be present. The file must exist and its filename must be a non-empty character string.

## Syntax

`[sys,inter]=gissmo2spinach(filename,subsystem)`

## Parameters / inputs

- `filename`: character string containing the GISSMO XML filename.
- `subsystem`: one-based index of the coupling matrix to import.

## Outputs

- `sys,inter`: Spinach data structures ready for `create`.

## Header notes

GISSMO supplies only chemical shifts, J-couplings, the non-selective linewidth, and the magnet field. Add further parameters by editing `sys` and `inter` as needed. As documented by `fwhm2rlx`, a linewidth-based relaxation rate is an approximation and should be treated as an upper bound when other broadening is present.
