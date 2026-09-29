# kernel/conventions/transforms/kelvin2hz.m

MATLAB source: [kernel/conventions/transforms/kelvin2hz.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/kelvin2hz.m)
Spinach Wiki: [kelvin2hz.m](https://spindynamics.org/wiki/index.php?title=kelvin2hz.m)

## Purpose

Converts temperature-valued energy scales, including Debye temperatures and thermal-energy scales in solid-state physics, to frequency values in hertz used in magnetic resonance.

## Usage

`hz = kelvin2hz(kelvin)`

## Input and output

- `kelvin`: a real numeric array of values in kelvin. Arrays of any dimensions are supported.
- `hz`: the corresponding frequency array in hertz, with the input array's shape.

The implementation checks that the input is numeric and real; it does not impose a particular array dimension.

## Conversion

The function applies `hz = 1.380649e-23 * kelvin / 6.62607015e-34`, i.e. `hz = (k_B / h) * kelvin`, using the exact SI values for the Boltzmann and Planck constants as written in the source. The conversion factor is approximately 20836619123.3 Hz/K. It returns ordinary frequency in hertz, not angular frequency.
