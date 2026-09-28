# kernel/pulses/vg_pulse.m

- Signature: `waveform=vg_pulse(pulse_name,npoints,duration)`

## Purpose

Generates Veshtort-Griffin shaped pulses from tabulated coefficients in the cited paper. Section 2.2 states that there are good reasons to believe these are the best possible pulses within their design specifications and basis sets.

## Numerical / algorithmic content

The function synthesizes cosine terms for `k=0:20` and sine terms for `k=1:20` on a grid spanning `0` to `2*pi`, then scales the waveform by `2*pi/duration`.

## Parameters / inputs

- `pulse_name` - character string selecting one of: `E0A`, `E0B`, `E100A`, `E100B`, `E200A`, `E200D`, `E200F`, `E300C`, `E300F`, `E400B`, `E300A`, `E500A`, `E500B`, `E500C`, `E600A`, `E600C`, `E600F`, `E800A`, `E800B`, or `E1000B`
- `npoints` - finite positive integer number of discrete time intervals in the pulse
- `duration` - finite positive real duration of the pulse, seconds

## Outputs

- `waveform` - pulse amplitude at each interval, with no phase modulation, normalized to produce a 90-degree pulse, rad/s

## Reference

https://doi.org/10.1002/cphc.200400018

Source Wiki page: https://spindynamics.org/wiki/index.php?title=vg_pulse.m
