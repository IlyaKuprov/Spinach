# Snapshot POD basis from supplied data

`pod_basis` makes an explicit weighted state-space basis for repeated simulations. Complex snapshot blocks and U*S factors are accepted in one common Euclidean coordinate system. Exact SVD minimises the supplied weighted Frobenius projection error; seeded randomised construction approximates that space and preserves the caller's RNG. Explicit numerical rank checks reject unsupported retained directions.

Weights multiply block energy, without automatic case normalisation. Diagnostics distinguish absolute/relative data residuals and exact versus sketch singular values. Pre-truncated factors define a surrogate data matrix; omitted original trajectories are not reconstructed. These diagnostics do not bound physical observable error or certify unseen parameter regions. Choose rank and training coverage explicitly, and use independent full-dynamics validation.
