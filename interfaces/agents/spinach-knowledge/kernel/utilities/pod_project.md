# Fixed-basis affine Liouvillian projection

`pod_project` prepares fixed affine Liouvillian terms, one initial state, and one detection vector in an explicitly supplied orthonormal Euclidean basis. All inputs must share state ordering and projectors. Its affine assembly identity is algebraic, not a statement of reduced-trajectory accuracy.

For Spinach, terms use L=H+1i*R+1i*K in rad/s; detection retains the Hermitian inner product. Non-Euclidean metrics require explicit whitening by the caller. No automatic basis choice, physical coordinate conversion, trace preservation, or positivity guarantee is implied. Preparation is distinct from repeated coefficient evaluation.

Liouvillian finiteness validation preserves sparse implicit zeros and scans dense terms without extracting a same-size double nonzero-entry vector.

Sparse bases and state vectors are screened without expanding their implicit zeros; sparse basis Gram matrices retain a sparse identity for orthogonality checking. Sparse coordinate-selection bases and sparse terms/states can therefore yield a sparse prepared model. Projected terms remain subject to computed-domain finiteness checks.
