# examples/cce/cce_polarized_bath.m

Signature: [answer,diagnostics]=cce_polarized_bath()

A two-proton illustrative bath with descending-m populations [0.8,0.2] is sampled using 64 Threefry states, seed 8800. Standard sys/inter/bas/parameters blocks call clexpn with hahn_echo, projected CCE order 2, and after-expansion averaging. The bath-state basis selector 'uniform' requests the standard Zeeman basis; the explicitly supplied probability rows are not uniform.

answer is the 21-by-1 normalised complex Hahn-echo signal on 0 to 200 microseconds. diagnostics retains the physical time axis, reduction, sampled configurations, estimator standard error, and arithmetic error bounds. The example plots the sampled signal only. Exhaustive product-density comparison and the original population array are preserved in test_cce_examples; no bath density matrices are assembled by the example.

This is not a material-geometry reproduction or a fitted thermal-relaxation model.

The related material calculation is Li et al. ([DOI](https://doi.org/10.1039/D6CP01591C), supplementary Fig. 7b). It uses a point electron at the pentacene substitution in naphthalene, a 40 Å proton bath, a 6 Å nuclear pair cutoff, and CCE-2 with single-proton density `(I+P*sigma_z)/2`. The toy example has the same product-population form at P=0.6, but it does not model the crystal, dipolar network, triplet transition, or the published echo curves. The public supplement does not specify the complete numerical input for regenerating those curves.
