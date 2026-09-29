# kernel/grids/sphtarea.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/sphtarea.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=sphtarea.m)

- Signature: `S=sphtarea(r1,r2,r3,sflag)`
- `r1`, `r2`, and `r3`: real three-element unit vectors for the triangle vertices. Each pairwise spherical arc length must be at most `pi/2`; each norm is checked against one within `sqrt(eps)`.
- `sflag`: optional character value `'signed'` or `'unsigned'`; omission defaults to `'unsigned'`.
- Output `S`: scalar area of the triangle on the unit sphere, numerically a solid angle in steradians. The signed result follows the vertex order; unsigned mode takes its absolute value. There are no time or frequency units.

After reshaping the vertices as columns, the source computes `2*atan2(det([r1 r2 r3]),dot(r1,r2)+dot(r2,r3)+dot(r3,r1)+1)`. This is the oriented spherical-triangle area rule; `'unsigned'` returns `abs(S)`. Inputs must be numeric and real with three elements, and the flag must match one of the two accepted character strings. The calculation is deterministic for fixed inputs. This geometry helper does not compute eigenfields, time evolution, or frequency offsets.
