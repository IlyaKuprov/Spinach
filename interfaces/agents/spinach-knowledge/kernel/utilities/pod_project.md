# Fixed-basis affine Galerkin projection

`pod_project(basis,terms,rho,coil)` prepares a fixed family in one state ordering: reduced terms are basis'*terms{k}*basis, initial coordinates are basis'*rho and detection coordinates are basis'*coil. Hermitian detection preserves coil'*rho. Liouvillians follow L=H+1i*R+1i*K in rad/s and propagation exp(-1i*L*t), with time in seconds.

Input data are power-of-two normalised in column slices bounded by reduced rank before each matrix product, then restored to physical units. Scaling does not materialise a second full dense Liouvillian. This avoids oversized intermediates when projected coordinates are finite. Nonzero components lost to normalisation and nonfinite final projected terms/states are unsupported and rejected. Finiteness checks preserve sparse implicit zeros; sparse bases use a sparse Gram identity. The dimension-dependent orthogonality budget is gamma(8*n)*r, where gamma(k)=k*eps/(1-k*eps), allowing the real/complex dot-product roundoff scale without repairing the supplied basis.

Projection commutes algebraically with affine assembly, but does not certify accurate reduced trajectories, positivity, trace preservation or stability. The caller supplies an orthonormal Euclidean basis, conformable finite double terms and state columns, common ordering and consistent projectors. No coordinate conversion or non-Euclidean metric is inferred. Training and evaluation remain separate. Fixed-basis offline/online projection follows Benner, Gugercin and Willcox, DOI:10.1137/130932715, section2.2.

Projection rejects operand ranges whose normalised real component products would underflow; it never silently substitutes a zero coordinate. This conservative range check does not certify componentwise conditioning.

The matrix-product range contract conservatively requires every nonzero operand-component pair to remain nonzero at the working scale; it can reject a physically finite result rather than silently lose a contribution.

The fixed Hermitian basis and its component minimum are prepared once for all projection slices, terms, and states. Range scans use bounded dense columns or stored sparse components; no full-basis range scan occurs inside the slice loop.

All-real projection operands retain real slice and prepared-array storage. Complex operands promote the result naturally; no imaginary component is discarded.

Physical restoration rejects loss of a nonzero real or imaginary working component, including partial loss in a still-nonzero complex coordinate. Counts before and after elementwise binary scaling identify that loss without retaining another projected array. Representable subnormal terms and states remain supported.

Normalised Galerkin column-slice sums retain TwoSum/cascaded corrections in projected-output-shaped workspace before physical restoration. A zero coordinate with lost correction information is rejected, rather than certifying a stationary generator. Compensation applies to computed slice products, not exact matrix multiplication. Independent dense/sparse and real/imaginary dyadic controls retain a projected 2^-56 term and its finite long-time phase.
