# kernel/orientation.m

- Signature: `H=orientation(Q,euler_angles)`

`orientation` evaluates the anisotropic Hamiltonian contribution at one orientation. For each rank represented in `Q`, it obtains a Wigner matrix from the three Euler angles and contracts its `(k,m)` coefficients with the corresponding operator components `Q{r}{k,m}`. Only nonzero components are accumulated. This produces the operator sum directly; it does not construct a separately rotated copy of the component tensor. The accumulator starts sparse, but the source does not convert the result back to sparse after the sum. The result is then Hermitian-symmetrised as `(H+H')/2`.

The input notes specify angles in radians relative to the input orientation and identify `Q` as the rotational basis returned by `hamiltonian.m`. The source check requires `Q` to be a cell array and the angles to be numeric, real, and contain three elements. It does not require a particular 1-by-3 shape or validate the nested dimensions of `Q` there. The source notes the same linear action can be used in Hilbert or Liouville space because the map `H -> [H, ]` is linear.

## References

- MATLAB source: [`kernel/orientation.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/orientation.m)
- Spinach Wiki: [`orientation.m`](https://spindynamics.org/wiki/index.php?title=orientation.m)
