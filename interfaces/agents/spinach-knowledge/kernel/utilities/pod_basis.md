# Snapshot POD basis from supplied data

`pod_basis` makes an explicit weighted state-space basis for repeated simulations. Complex snapshot blocks and U*S factors are accepted in one common Euclidean coordinate system. Exact SVD minimises the supplied weighted Frobenius projection error; seeded randomised construction approximates that space and preserves the caller's RNG. Explicit numerical rank checks reject unsupported retained directions.

Weights multiply block energy, without automatic case normalisation. Diagnostics distinguish absolute/relative data residuals and exact versus sketch singular values. Pre-truncated factors define a surrogate data matrix; omitted original trajectories are not reconstructed. These diagnostics do not bound physical observable error or certify unseen parameter regions. Choose rank and training coverage explicitly, and use independent full-dynamics validation.

Randomised construction applies weighted sparse or dense blocks directly, without forming a dense concatenated snapshot matrix. Range, adjoint-range, and compressed-data storage scale with the sketch rank; energy and residual measurement use column slices no wider than the sketch rank. Both methods omit unused right singular vectors.
