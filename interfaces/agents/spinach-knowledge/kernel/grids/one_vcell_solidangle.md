# kernel/grids/one_vcell_solidangle.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/one_vcell_solidangle.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=one_vcell_solidangle.m) · [method reference](https://doi.org/10.1109/TBME.1983.325207)

- Signature: `S=one_vcell_solidangle(v,centre)`
- `v`: finite real `3 x n` matrix; each column is a unit-vector polygon vertex, with squared norm within `1e-6` of one.
- `centre`: optional finite real `3 x 1` unit vector, with norm within `1e-6` of one.
- Output `S`: scalar signed sum for the convex spherical polygon, reported by the source in radians; geometrically it is unit-sphere solid angle (steradians). The sign follows vertex ordering. No time or frequency units occur.

Without `centre`, the source makes a fan of triangles `[v(:,1),v(:,k),v(:,k+1)]` for `k=2,...,n-1`. With `centre`, it closes the sequence by appending the first vertex and makes one triangle `[centre,v(:,k),v(:,k+1)]` per edge. For each triangle matrix `T`, the contribution is `2*atan2(det(T),1+T1·T2+T2·T3+T3·T1)`; `S` is the sum of these oriented contributions. Thus the centre changes the fan triangulation, not the output unit. The source checks numeric type, reality, finiteness, row count, and unit lengths; it does not test convexity or enforce a minimum vertex count. The calculation is deterministic for fixed inputs. This geometry helper does not compute eigenfields, time evolution, or frequency offsets.
