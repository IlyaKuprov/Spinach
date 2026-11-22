# kernel/utilities/conmat.m

The call `conmat(xyz,r0)` returns symmetric sparse logical connectivity with an empty diagonal. Each Cartesian separation must be strictly less than the cutoff and the Euclidean length must not exceed it. Distinct coincident points are connected. Occupied spatial bins limit the candidate search without allocating a dense pair mask.

## References

- Source: [kernel/utilities/conmat.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/conmat.m)
- Spinach Wiki: [conmat.m](https://spindynamics.org/wiki/index.php?title=conmat.m)
