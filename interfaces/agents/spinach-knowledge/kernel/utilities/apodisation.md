# kernel/utilities/apodisation.m

- Signature: `fid=apodisation(spin_system,fid,winfuns,fp_half)`

## Purpose

Apply independent window functions to the active dimensions of a free induction decay (FID). By default, the first point along each active, non-singleton dimension is halved to meet Fourier-transform symmetry requirements.

## Parameters / inputs

- `spin_system` - Spinach spin-system structure used for reporting.
- `fid` - numeric FID array. The time origin is at the first point in each dimension (for example, index `(1,1)` in 2D and `(1,1,1)` in 3D).
- `winfuns` - cell array in the form `{{spec},{spec},...}`, with one specification for each non-singleton dimension in dimension order; singleton dimensions are omitted. An empty specification `{}` leaves that dimension inactive. Supported specifications:
  - `{'none'}` - unit window; the first point is still halved when `fp_half` is true.
  - `{'crisp'}` - `cos(x)^8`, with `x` from 0 to `pi/2`.
  - `{'exp',k}` - `exp(-k*x)`, with `x` from 0 to 1.
  - `{'gauss',k}` - `exp(-k*x^2)`, with `x` from 0 to 1.
  - `{'cos'}` - `cos(x)`, with `x` from 0 to `pi/2`.
  - `{'sin'}` - `sin(x)`, with `x` from 0 to `pi`.
  - `{'sqcos'}` - `cos(x)^2`, with `x` from 0 to `pi/2`.
  - `{'sqsin'}` - `sin(x)^2`, with `x` from 0 to `pi`.
  - `{'kaiser',k}` - Kaiser window with side-lobe attenuation factor `k` and a peak at the centre of the FID.
  - `{'bad-z1',k}` - emulation of a misset Z1 shim; `k` is dimensionless and proportional to shim current. The source gives 10 as a good guess for 1H NMR at 600 MHz.
  - `{'bad-z2',k}` - emulation of a misset Z2 shim; `k` is dimensionless and proportional to shim current. The source gives 40 as a good guess for 1H NMR at 600 MHz.
- `fp_half` - optional logical flag; defaults to true. Set false to skip first-point halving, for example when applying successive windows to the same dimension.

## Numerical content

The function first applies the optional first-point factors, then constructs and multiplies each requested one-dimensional window along its corresponding array dimension. The shim-emulation windows use the source's sinc and Fresnel-integral constructions, respectively. A window parameter `k` must be a finite real scalar.

## Output

- `fid` - the apodised FID array.