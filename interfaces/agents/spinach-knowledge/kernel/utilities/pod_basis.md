# Snapshot POD basis from supplied data

`pod_basis` makes an explicit weighted state-space basis for repeated simulations. Complex snapshot blocks and U*S factors are accepted in one common Euclidean coordinate system. Exact SVD minimises the supplied weighted Frobenius projection error; seeded randomised construction approximates that space and preserves the caller's RNG. Explicit numerical rank checks reject unsupported retained directions.

Weights multiply block energy, without automatic case normalisation. Diagnostics distinguish absolute/relative data residuals and exact versus sketch singular values. Pre-truncated factors define a surrogate data matrix; omitted original trajectories are not reconstructed. These diagnostics do not bound physical observable error or certify unseen parameter regions. Choose rank and training coverage explicitly, and use independent full-dynamics validation.

Zero data is detected from the unsquared weighted norm. Squared `total_energy` can round to zero for tiny nonzero snapshots; `residual_rel` divides by the unsquared norm. Energy overflow remains an error, and this handling does not change numerical rank tolerances or certify arbitrary subnormal inputs.

Randomised construction applies weighted sparse or dense blocks directly, without forming a dense concatenated snapshot matrix. For n state coordinates, M positively weighted columns, and working sketch rank s, it retains n-by-s range arrays and M-by-s or s-by-M sketch/compressed arrays; total workspace therefore also depends on M. Energy and residual measurement use column slices no wider than the sketch rank. Both methods omit unused right singular vectors and request singular values directly as a vector, without allocating a diagonal singular-value matrix.

Zero-weight blocks remain subject to the input contract but do not enter construction or residual arithmetic. Internal sketches cannot exceed the positively weighted column count. Returned singular-value vectors retain their original supplied-data lengths with zero padding; numerical rank tolerance retains the original supplied-data dimensions.

Randomised sketches and returned bases remain real when every positively weighted block is real; a zero-weight complex block does not force imaginary storage. Active complex data retain complex Gaussian sketches. Both cases preserve caller RNG and fixed-seed repeatability.

Sketch coefficients are released before compressed snapshots are allocated. After numerical rank validation, randomised mode construction multiplies only the retained compressed left singular vectors; compressed data and unused range/mode arrays are released before residual measurement. Returned singular-value diagnostics are unchanged.
