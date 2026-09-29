# etc/estimators/guess_csa_pro.m

- MATLAB implementation: [etc/estimators/guess_csa_pro.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/estimators/guess_csa_pro.m)

`CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)`

This auxiliary function, called by the `protein.m` import module (direct calls are discouraged), assigns approximate CSA tensors to amide `15N`, carbonyl `13C`, and amide H atoms from local protein geometry. The output `CSAs` is a cell array with one slot per input atom; assigned entries are 3×3 tensors in ppm, and atoms without a supported assignment remain empty. These are rough guesses, not relaxation-quality tensors; use experimentally or otherwise justified tensors for accurate relaxation analysis.

A direct call needs four inputs. Supply one aligned entry per atom in `aa_nums`, `pdb_ids`, and `coords`; the latter two are cell arrays of PDB atom labels and coordinate vectors. `options` must be a structure; an empty structure selects the default:

```matlab
options=struct();                 % defaults to options.nh_csa='tcb'
CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options);
```

`options.nh_csa` controls the eigenvalues for both amide N and amide H: `'tcb'` (default), `'bax'`, or `'pol'`. Values are listed in XX, YY, ZZ order (ppm); the carbonyl values do not depend on this option.

| Site | `'tcb'` | `'bax'` | `'pol'` |
|---|---:|---:|---:|
| Amide 15N | −125, 45, 80 | −108, 62, 46 | −92.4, 34.7, 57.7 |
| Amide H | 7, 0, −7 | 6, 0, −6 | 6.66, 0.67, −7.33 |
| Carbonyl 13C | 70, 5, −75 | 70, 5, −75 | 70, 5, −75 |

## Geometry-to-tensor mapping

For each supported site the tensor is assembled as `V*D*V'`, where `D` is the chosen principal-value diagonal and `V` contains axes built from normalised bond directions and cross products. For the amide N in residue `n+1`, the routine uses the preceding residue's carbonyl C together with the current N and H: ZZ is collinear with the C–N bond, YY is normal to the C–N–H plane, and XX completes the frame. The carbonyl C in residue `n` requires that residue's C and CA plus N in residue `n+1`; XX follows C-to-CA and ZZ is normal to the carbonyl plane. Amide H in residue `n` requires H, N, and CA in that residue and C in residue `n−1`; its YY axis follows H-to-CA, XX is normal to the peptide plane, and ZZ completes the frame.

This adjacency encodes the assumption that amino-acid numbers progress from N- to C-terminus. The source explicitly says proline is not handled. Missing required atoms cause that site to be skipped (with a message); distance checks use thresholds of 2.0 for the amide-N and carbonyl constructions, and 1.2 (N–H), 1.6 (N–CA), and 1.5 (N–preceding-C) for amide H. A failed distance check raises an “Amino acid numbering is not sequential” error. The source does not state coordinate units, so these thresholds apply in the same length units as the supplied coordinates; the input checks do not validate vector contents or units.

## Sources

The TCB choices cite [doi:10.1007/s10858-006-9037-6](https://doi.org/10.1007/s10858-006-9037-6), [doi:10.1021/ja00083a028](https://doi.org/10.1021/ja00083a028), and [doi:10.1021/ja042863o](https://doi.org/10.1021/ja042863o). The Bax branches cite [doi:10.1021/ja0016194](https://doi.org/10.1021/ja0016194); the Case–Polenova–Gronenborn branches cite [doi:10.1039/C8CP00647D](https://doi.org/10.1039/C8CP00647D).

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=Guess_csa_pro.m).
