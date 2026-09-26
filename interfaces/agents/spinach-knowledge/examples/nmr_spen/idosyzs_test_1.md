# examples/nmr_spen/idosyzs_test_1.m

- Signature: `idosyzs_test_1()`

## Purpose

Diffusion attenuation during soft pulses in a simplified model sequence of the Zangger-Sterk pure shift iDOSY with fitting using a modified ver- sion of the Stejskal Tanner equation, as described in: Calculation time: seconds on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Diffusion attenuation during soft pulses in a simplified model sequence
- of the Zangger-Sterk pure shift iDOSY with fitting using a modified ver-
- sion of the Stejskal Tanner equation, as described in:
- Calculation time: seconds on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Isotopes
- Chemical shift
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Diffusion coefficient
