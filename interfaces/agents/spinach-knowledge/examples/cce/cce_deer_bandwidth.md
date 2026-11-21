# examples/cce/cce_deer_bandwidth.m

Signature: `[answer,diagnostics]=cce_deer_bandwidth()`

The self-contained three-electron toy system uses ordinary sys/inter input, no full-bath basis, and one call to clexpn with @bath_pumped_hahn. It is a finite-cluster demonstration, not a published material fit.

The observed spin is physical site 1, with transition ms=[0.5 -0.5]. The pump selects bath site 2 at zero frequency, with 100 kHz halfwidth (200 kHz full width), Hz units, fraction 0.25, and allow_empty=true. The ideal observed-spin refocusing pulse is at half the total echo duration. The projected expansion retains both bath spins at order two. Numerical tolerances remain atol=rtol=1e-8 and the memory budget is 2^28 bytes.

The default no-mean-field product-density ensemble omits a selected exterior spin from each retained marginal cluster, including the empty cluster. Static conditioned exterior pumping is explicitly unsupported. Changing the spectral window to select no spins removes the bath rotation when allow_empty=true.

The first output is the natural complex echo array, and the second contains native reduction, ensemble, and raw numerical diagnostics. test_cce_examples checks the result against direct three-spin density-matrix propagation; test_clexpn_context checks the analytic one-bath-spin response and pump selector validation.
