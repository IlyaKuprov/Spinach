# kernel/grids/sphtarea.m

- Signature: `S=sphtarea(r1,r2,r3,sflag)`

## Purpose

Compute the area of a curvilinear triangle on the unit sphere from its three vertex coordinates.

## Inputs

- `r1`, `r2`, `r3`: three-element real unit vectors giving the Cartesian coordinates of the vertices. Each pairwise arc length must not exceed `pi/2`.
- `sflag`: `'signed'` accounts for surface-normal direction; `'unsigned'` returns a nonnegative area and is the default.

## Output

- `S`: spherical triangle surface area.

## Method

The function reshapes the vertices into columns and computes `S=2*atan2(det([r1 r2 r3]),dot(r1,r2)+dot(r2,r3)+dot(r3,r1)+1)`. For `'unsigned'`, it returns `abs(S)`.

Source link: <https://spindynamics.org/wiki/index.php?title=sphtarea.m>