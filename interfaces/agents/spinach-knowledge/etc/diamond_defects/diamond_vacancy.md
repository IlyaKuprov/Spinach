# etc/diamond_defects/diamond_vacancy.m

`[sys,inter]=diamond_vacancy(parameters)`

Builds a one-electron Spinach model for the requested vacancy-family centre. The implementation supplies the centre-specific g and zero-field-splitting (ZFS) principal values below; ZFS values are in MHz. `w29` uses electron label `E4`; all other listed centres use `E3`.

| Centre | g principal values | ZFS principal values (MHz) |
|---|---|---|
| `r4_w6`, `w6`, `r4` | 2.0022, 2.0026, 2.0013 | 105, 197, −303 |
| `w29` | 2.002, 1.997, 2.005 | 297, 156, −453 |
| `r5` | 2.00275, 2.00265, 2.00205 | 283, 244, −524 |
| `o1` | 2.00299, 2.00273, 2.00212 | 109, 95, −205 |
| `r6` | 2.00299, 2.00273, 2.00212 | 62, 59, −120 |
| `r10` | 2.00295, 2.00269, 2.00212 | 36, 36, −73 |
| `r11` | 2.00301, 2.00278, 2.00200 | 27, 27, −53 |

## Parameters

- `parameters.centre`: `'r4_w6'`, `'w6'`, `'r4'`, `'w29'`, `'r5'`, `'o1'`, `'r6'`, `'r10'`, or `'r11'`.
- `parameters.orientation`: `'111'`, `'110'`, or `'100'`; the corresponding crystal plane normal is aligned with the magnetic-field axis.

The routine rotates the selected tensors into the requested orientation, converts the ZFS tensor to Spinach's interaction representation, and returns the system and interaction structures. Both parameters must be fields of a structure, and the centre and orientation must be character values. The implementation reports an error for an unrecognised centre or orientation.

## Sources

- R4/W6: Twitchen et al., *Physical Review B* **59**, 12900 (1999), [doi:10.1103/PhysRevB.59.12900](https://doi.org/10.1103/PhysRevB.59.12900).
- W29: Kirui et al., *Diamond and Related Materials* **8**, 1569 (1999), [doi:10.1016/S0925-9635(99)00037-0](https://doi.org/10.1016/S0925-9635(99)00037-0).
- R5/O1/R6/R10/R11: Iakoubovskii and Stesmans, *Physical Review B* **66**, 045406 (2002), [doi:10.1103/PhysRevB.66.045406](https://doi.org/10.1103/PhysRevB.66.045406). The ZFS table values are cross-checked against Ball, PhD thesis, OIST Graduate University (2021).

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_vacancy.m).
