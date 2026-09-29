# etc/diamond_defects/diamond_vacancy.m

- MATLAB implementation: [etc/diamond_defects/diamond_vacancy.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_vacancy.m)

`[sys,inter]=diamond_vacancy(parameters)`

Build a one-electron Spinach model for one of the supported diamond vacancy-family centres. For example:

```matlab
parameters.centre='r4';
parameters.orientation='111';
[sys,inter]=diamond_vacancy(parameters);
```

The principal values below are the source's g values (dimensionless) and zero-field-splitting (ZFS) values in MHz, in the listed principal-axis order.

| Centre label(s) | g principal values | ZFS principal values (MHz) |
|---|---|---|
| `r4_w6`, `w6`, `r4` | 2.0022, 2.0026, 2.0013 | 105, 197, −303 |
| `w29` | 2.002, 1.997, 2.005 | 297, 156, −453 |
| `r5` | 2.00275, 2.00265, 2.00205 | 283, 244, −524 |
| `o1` | 2.00299, 2.00273, 2.00212 | 109, 95, −205 |
| `r6` | 2.00299, 2.00273, 2.00212 | 62, 59, −120 |
| `r10` | 2.00295, 2.00269, 2.00212 | 36, 36, −73 |
| `r11` | 2.00301, 2.00278, 2.00200 | 27, 27, −53 |

## Inputs and constraints

`parameters` must be a structure containing character-vector fields `.centre` and `.orientation`. The accepted centres are exactly the table labels; the centre is lowercased before matching. The accepted orientations are `'111'`, `'110'`, and `'100'`. Each names a crystal-plane normal to align with the magnetic-field axis; it is not a free set of Euler angles.

## How the model is oriented

The routine constructs centre-specific principal-axis frames (separate frames for R4/W6, W29, and the vacancy-chain centres), forms the g and ZFS matrices from the tabulated diagonal principal values, then rotates them for the requested crystal orientation. It assigns electron isotope `E4` only for `w29`; every other supported centre uses `E3`. The ZFS values are scaled from MHz to Hz in the source and converted with `mat2ias` into Spinach's interaction representation before being placed in `inter.coupling.matrix{1,1}`. The rotated g matrix is returned in `inter.zeeman.matrix{1}`.

The outputs are Spinach system and interaction structures for this single electron. The function does not add surrounding nuclei or other couplings; those must be supplied separately if needed.

## Parameter sources

R4/W6 values are attributed to Twitchen et al., *Physical Review B* **59**, 12900 (1999), [doi:10.1103/PhysRevB.59.12900](https://doi.org/10.1103/PhysRevB.59.12900). W29 values are attributed to Kirui et al., *Diamond and Related Materials* **8**, 1569 (1999), [doi:10.1016/S0925-9635(99)00037-0](https://doi.org/10.1016/S0925-9635(99)00037-0). R5/O1/R6/R10/R11 values are attributed to Iakoubovskii and Stesmans, *Physical Review B* **66**, 045406 (2002), [doi:10.1103/PhysRevB.66.045406](https://doi.org/10.1103/PhysRevB.66.045406). The source says the ZFS table was cross-checked against Ball's 2021 OIST PhD thesis; that attribution is retained, not presented as an independent validation here.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_vacancy.m).
