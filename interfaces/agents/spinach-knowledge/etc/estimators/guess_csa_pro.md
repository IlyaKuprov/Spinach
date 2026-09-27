# etc/estimators/guess_csa_pro.m

`CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)`

This auxiliary routine estimates approximate chemical-shift-anisotropy (CSA) tensors in ppm for amide ¹⁵N, carbonyl ¹³C, and amide H atoms from a protein's local coordinates. It is called by the `protein.m` import module; direct calls are discouraged. Residue numbering is assumed to run from the N-terminus to the C-terminus, and proline is not handled. Missing required atoms leave the corresponding output cell empty; successful assignments and failures are reported to the command window. The estimates are explicitly described as very approximate; accurate relaxation analysis requires user-supplied tensors.

## Inputs and output

- `aa_nums`: amino-acid number for each input atom.
- `pdb_ids`: cell array of PDB atom identifiers; `coords`: matching cell array of coordinate vectors. The three input arrays must have the same number of elements.
- `options.nh_csa`: selects the nitrogen and proton eigenvalues: `'tcb'` (default), `'bax'`, or `'pol'`. The carbonyl-carbon values do not depend on this option.
- `CSAs`: cell array in input order; assigned elements contain a 3-by-3 tensor, while unassigned elements remain empty. (The source comment names the output `CSAa`, but the function signature and returned variable are `CSAs`.)

Each tensor is formed as `V*D*V'`, with the geometry-derived principal axes in `V` and the eigenvalues in `D`. Eigenvalues below are listed in XX, YY, ZZ order, in ppm.

| Site | `'tcb'` | `'bax'` | `'pol'` |
|---|---:|---:|---:|
| Amide ¹⁵N | −125, 45, 80 | −108, 62, 46 | −92.4, 34.7, 57.7 |
| Amide H | 7, 0, −7 | 6, 0, −6 | 6.66, 0.67, −7.33 |

The carbonyl ¹³C values are 70, 5, −75 ppm (XX, YY, ZZ) for all three options.

## Geometrical assignments

- The nitrogen in residue `n+1` is assigned from the preceding residue's carbonyl C and its own N and H. The ZZ axis follows the N-to-preceding-C direction; YY is normal to the plane defined with N-to-H, and XX completes the frame.
- The carbonyl C in residue `n` requires C and CA in residue `n` and N in residue `n+1`. XX follows C-to-CA; ZZ is normal to the plane including the next N, and YY completes the frame.
- Amide H in residue `n` requires H, N, and CA in `n`, plus C in `n−1`. YY follows CA-to-H; XX is normal to the plane defined by CA-to-N and preceding-C-to-N, and ZZ completes the frame.

The routine checks required input container types and matching element counts. It also uses geometry-distance checks to flag apparent non-sequential residue numbering; the source does not state the coordinate distance units.

The tensor orientations are associated in the source with [doi:10.1007/s10858-006-9037-6](https://doi.org/10.1007/s10858-006-9037-6), [doi:10.1021/ja00083a028](https://doi.org/10.1021/ja00083a028), and [doi:10.1021/ja042863o](https://doi.org/10.1021/ja042863o). The `'bax'` and `'pol'` branches also cite [doi:10.1021/ja0016194](https://doi.org/10.1021/ja0016194) and [doi:10.1039/C8CP00647D](https://doi.org/10.1039/C8CP00647D), respectively.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=Guess_csa_pro.m).
