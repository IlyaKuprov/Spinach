# kernel/line_shapes/lorentzfun.m

- Signature: `[real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)`

## Purpose

Compute a Lorentzian line shape in magnetic-resonance notation with phase distortion. The absorption mode integrates to `ampl/2` under the one-sided FID Fourier-transform convention; the source notes that `lorentzcon()` integrates to `ampl`.

## Physical / mathematical content

Set `gam=fwhm/2`, `u=(x-offs)/gam`, and `L=ampl/(2*pi*gam)/(1+u^2)`. The outputs are `real_part=L*cos(phi)-u*L*sin(phi)` and `imag_part=L*sin(phi)+u*L*cos(phi)`. At `phi=0`, the real component is the absorptive Lorentzian and the imaginary component is the dispersive component. The width parameter `fwhm` is the full width at half-maximum.

## Numerical / algorithmic content

The function validates the inputs, sets `gam=fwhm/2`, and evaluates the phase-mixed components elementwise over `x`. Both returned arrays have the same size as `x`.

## Parameters / inputs

- `offs` - real scalar peak offset from zero
- `ampl` - real scalar amplitude multiplier
- `fwhm` - positive real scalar full width at half-maximum
- `x` - numeric real array of any dimension
- `phi` - real scalar phase distortion in radians

## Outputs

- `real_part` - real component of the phase-distorted line shape, same size as `x`
- `imag_part` - imaginary component of the phase-distorted line shape, same size as `x`

## Implementation structure

A private consistency check rejects nonnumeric or nonreal `x`, a nonpositive or nonscalar `fwhm`, and nonscalar or nonnumeric/nonreal `offs`, `ampl`, or `phi`. The calculation then uses `gam=fwhm/2` in the Lorentzian denominator and applies the phase rotation to the two components.

## Reference

- [Spin Dynamics documentation for `lorentzfun.m`](https://spindynamics.org/wiki/index.php?title=lorentzfun.m)
