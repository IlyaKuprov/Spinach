# kernel/conventions/transforms/hartree2joule.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/hartree2joule.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hartree2joule.m)

## Signature

`energy=hartree2joule(energy)`

## Purpose and conversion

Multiplies an energy value expressed in Hartree by the source's fixed conversion factor:

`energy_out=2625499.62*energy_in`

The source describes a Hartree as twice the ground-state ionisation energy of hydrogen. The factor maps one Hartree to `2625499.62` J/mol. Although the output-parameter comment says Joules, the implemented conversion is molar energy, not joules per particle.

## Inputs and output

- `energy`: real numeric array; input unit is Hartree.
- Output `energy`: converted numeric array in J/mol, with the same array layout for supported inputs.

The implementation rejects nonnumeric or nonreal input. It imposes no explicit size or dimensionality check.

## Reference

- [Spinach Wiki: hartree2joule.m](https://spindynamics.org/wiki/index.php?title=hartree2joule.m)
