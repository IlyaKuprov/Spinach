# etc/diamond_defects/diamond_nvm_gs.m

- Signature: [sys,inter]=diamond_nvm_gs(parameters)

## Purpose

Constructs the NV centre ground-state spin system in diamond using the magnetic parameters of Felton et al., *Phys. Rev. B* **79**, 075203 (2009), https://doi.org/10.1103/PhysRevB.79.075203.

## Physical / mathematical content

The electron is represented as E3. Its principal g values are 2.0031, 2.0031, and 2.0029, with axial zero-field splitting D = 2872 MHz. For 14N, the electron–nuclear hyperfine principal values are -2.70, -2.70, and -2.14 MHz, and the nitrogen quadrupole interaction has D = -5.01 MHz. For 15N, the routine uses hyperfine values 3.65, 3.65, and 3.03 MHz and adds no quadrupole interaction.

## Numerical / algorithmic content

The tensors are defined in the trigonal principal-axis frame and rotated to align the selected crystal direction with the field. If omitted, orientation and nitrogen isotope default to '111' and '14N', respectively.

## Parameters / inputs

- parameters.orientation: '111', '110', or '100'; default '111'.
- parameters.nitrogen: '14N' or '15N'; default '14N'.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

The routine fills the isotope list and Zeeman/coupling matrices, including the electron zero-field splitting and the isotope-dependent nitrogen interactions. The consistency check validates the structure and, when supplied, the orientation's character-string type; the orientation switch rejects unsupported values.
