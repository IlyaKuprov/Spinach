# etc/molecules/cyprinol.m

- MATLAB implementation: [etc/molecules/cyprinol.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/cyprinol.m)

## Call and returned model

`[sys,inter,bas] = cyprinol()` takes no inputs. It constructs a fixed 69-spin cyprinol model: 42 `1H` spins followed by 27 `13C` spins. It returns the isotope list and chemical shifts in `sys` and `inter`, a scalar-coupling array in `inter`, and methyl permutation-symmetry settings in `bas`. The source does not assign a coordinate array or proton–proton J-couplings (both are marked TODO).

## Parameter provenance and model limits

The source cites [10.1002/mrc.4782](https://doi.org/10.1002/mrc.4782) for isotropic chemical shifts and J couplings, then explicitly says values absent from that report were estimated by “tossing a twenty-sided coin.” Treat those missing-data values as placeholders, not validated physical estimates. The implementation uses repeated heuristic coupling constants (including 40, 3, and 0.3 in the C–C block and 150 and 2 in the C–H block); the source does not state their units or provide a calibration procedure. Numeric shifts are included in the source, but it does not explain how individual assignments map to the cited paper.

Three proton triplets are declared as `S3` symmetry sets in `bas.sym_spins`: `H21a/H21b/H21c`, `H18a/H18b/H18c`, and `H19a/H19b/H19c`. These are the model's explicit permutation-symmetry instructions; they are not a claim that all sites in the molecule are equivalent under every experimental condition. For a richer test spin system, the source recommends `strychnine.m`.

Source: [cyprinol.m in Spinach](https://spindynamics.org/wiki/index.php?title=cyprinol.m); cited paper DOI: [10.1002/mrc.4782](https://doi.org/10.1002/mrc.4782).
