# kernel/pulses/vg_pulse.m

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/vg_pulse.m
- Wiki: https://spindynamics.org/wiki/index.php?title=vg_pulse.m
- Signature: `waveform=vg_pulse(pulse_name,npoints,duration)`

## Purpose

Generates the tabulated Veshtort–Griffin shaped pulses associated with the cited paper. The source repeats the paper's statement that there are good reasons to believe these are the best possible pulses within their design specifications and basis sets (Section 2.2); that is the paper's claim, not a general guarantee for other designs.

## Coefficients and sampling

The available pulse names select columns of coefficient tables embedded in the MATLAB source; the function does not take a waveform-file path or read an external pulse file. For the selected name, it sums cosine harmonics `k=0:20` and sine harmonics `k=1:20` on a `npoints`-element grid from `0` through `2*pi`, then applies the factor `2*pi/duration`. The output is sampled amplitude, not a phase-modulated waveform.

## Inputs and output

- `pulse_name` - character string selecting one of `E0A`, `E0B`, `E100A`, `E100B`, `E200A`, `E200D`, `E200F`, `E300C`, `E300F`, `E400B`, `E300A`, `E500A`, `E500B`, `E500C`, `E600A`, `E600C`, `E600F`, `E800A`, `E800B`, or `E1000B`
- `npoints` - finite positive integer number of discrete pulse intervals
- `duration` - finite positive real pulse duration, in seconds
- `waveform` - amplitude at the sampled intervals, in radians per second; the source documentation describes it as having no phase modulation and being normalised to produce a 90-degree pulse

There is no filter or phase input: pulse identity, point count, and duration are the controls exposed by this function.

## Reference

- Veshtort–Griffin pulse design paper: https://doi.org/10.1002/cphc.200400018
