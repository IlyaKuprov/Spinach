# Fixed-basis affine Liouvillian projection

`pod_project` prepares fixed affine Liouvillian terms, one initial state, and one detection vector in an explicitly supplied orthonormal Euclidean basis. All inputs must share state ordering and projectors. Its affine assembly identity is algebraic, not a statement of reduced-trajectory accuracy.

For Spinach, terms use L=H+1i*R+1i*K in rad/s; detection retains the Hermitian inner product. Non-Euclidean metrics require explicit whitening by the caller. No automatic basis choice, physical coordinate conversion, trace preservation, or positivity guarantee is implied. Preparation is distinct from repeated coefficient evaluation.

Liouvillian finiteness validation preserves sparse implicit zeros and scans dense terms without extracting a same-size double nonzero-entry vector.
