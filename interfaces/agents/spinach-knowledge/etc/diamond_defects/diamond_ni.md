# etc/diamond_defects/diamond_ni.m

- Signature: [sys,inter]=diamond_ni(parameters)

## Purpose

Constructs a spin system for the listed nickel-related diamond centres. The W8 magnetic parameters are from Isoya, Kanda, Norris, Tang, and Bowman, *Phys. Rev. B* **41**, 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905. Parameters for the other nickel centres are from Nadolinny et al., *Crystals* **7**, 237 (2017), https://doi.org/10.3390/cryst7080237.

The W8 quartet entry uses zero zero-field splitting because the cited parameter set reports none; off-central transitions are treated as unresolved, not explicitly modelled. In the Nadolinny review, the Table 2 zero-field-splitting column is in tesla, as its heading states. The source also records reported values of D = -171 GHz for NOL1/NIRIM5 and D = 31.72 GHz for AB5; both exceed the 9.5 GHz X-band quantum, while the AB5 splitting is below the 94 GHz W-band quantum used in the shipped example.

## Physical / mathematical content

Depending on the selected centre, the system contains an electron (or an effective E3/E4 spin), optional nickel and carbon nuclei, and the tabulated nitrogen hyperfine tensors. The routine includes the centre's electron g tensor and, where specified, its zero-field-splitting tensor. For W8, the nickel isotope can be omitted with 'none'; 61Ni is assigned its tabulated hyperfine coupling.

## Numerical / algorithmic content

Tabulated principal values are converted to frequency units where necessary and transformed from the centre frame to the selected crystal orientation. For W8, up to four nearest-neighbour 13C nuclei can be included. They are placed in source order on [111], [1-1-1], [-11-1], and [-1-11] bonds. Selecting fewer than four therefore selects a particular isotopomer, not an average: at orientation '111', the [111] carbon splits by 1.339 mT and each other listed carbon by 0.451 mT. To represent a mixture, include all four and select or weight isotopomers in the calling script.

## Parameters / inputs

- parameters.centre: 'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5', 'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5', 'nol1', or 'nirim5'.
- parameters.orientation: '111', '110', or '100'; the corresponding crystal-plane normal is aligned with the magnetic field.
- parameters.nickel: required for W8; '61Ni', 'none', or another isotope string. It is ignored for other centres.
- parameters.n_13c: required for W8; integer from 0 to 4. It is only supported for W8 and must be zero otherwise.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

After validating the input, the function selects centre-specific tensors and nuclei, constructs the orientation rotation, then fills the Zeeman and coupling matrices. The W8 path optionally appends the requested nearest-neighbour 13C nuclei; the NOL1/NIRIM5 and AB5 paths include their zero-field-splitting tensors.
