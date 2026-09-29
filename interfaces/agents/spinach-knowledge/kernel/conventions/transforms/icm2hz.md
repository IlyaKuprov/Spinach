# kernel/conventions/transforms/icm2hz.m

MATLAB source: [kernel/conventions/transforms/icm2hz.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/icm2hz.m)
Spinach Wiki: [icm2hz.m](https://spindynamics.org/wiki/index.php?title=icm2hz.m)

## Purpose

Converts inverse-centimetre values used in spectroscopy to frequency values in hertz, the units preferred in magnetic resonance.

## Usage

`hz = icm2hz(icm)`

## Input and output

- `icm`: a real numeric array of values in inverse centimetres. Arrays of any dimensions are supported.
- `hz`: the corresponding array in hertz, with the input array's shape.

The implementation checks that the input is numeric and real; it does not impose a particular array dimension.

## Conversion

The function applies `hz = 100 * 299792458 * icm`. The factor 100 converts inverse centimetres to inverse metres, and 299792458 m/s is the exact speed of light used by the source. Thus 1 cm^-1 is converted to 29979245800 Hz. This is a frequency conversion, not an angular-frequency conversion.
