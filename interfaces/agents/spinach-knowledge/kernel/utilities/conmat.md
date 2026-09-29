# kernel/utilities/conmat.m

## Purpose

Computes the molecular connectivity matrix for a set of particles given their Cartesian coordinates and a connection distance cutoff. The header comment states that the algorithm has N·log(N) asymptotic complexity scaling with respect to the number of atoms.

## Behaviour

- Syntax: `conmatrix=conmat(xyz,r0)`.
- Input consistency is enforced by a local `grumble` subfunction, which errors if `xyz` is not a real numeric matrix with three columns, or if `r0` is not a positive real number.
- The function sorts the particles by each Cartesian coordinate in turn (X, then Y, then Z). For each sorted axis, it scans pairs and marks a connection in a logical matrix (`A` for X, `B` for Y, `C` for Z) when the absolute difference of the sorted coordinates along that axis is less than `r0`; the inner loop breaks at the first pair whose coordinate difference reaches `r0`.
- Each axis matrix is symmetric: both `(x_index(n),x_index(k))` and `(x_index(k),x_index(n))` positions are set (and analogously for Y and Z).
- The three axis matrices are combined with logical AND (`conmatrix=A&B&C`) to produce a candidate ("dirty") matrix.
- A cleanup pass then examines each candidate pair found via `find` and removes the entry if the Euclidean distance `norm(xyz(row(n),:)-xyz(col(n),:),2)` exceeds `r0`.
- The final result is converted with `sparse` before being returned.

## Inputs and outputs

**Inputs**

- `xyz` — an array with N rows and three columns giving the Cartesian coordinates of each particle; must be real and numeric.
- `r0` — the distance below which particles are considered "connected"; must be a positive real number.

**Output**

- `conmatrix` — a sparse logical matrix containing 1 at positions corresponding to connected particle pairs.

## References

- Source: [kernel/utilities/conmat.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/conmat.m)
- Spinach Wiki: [conmat.m](https://spindynamics.org/wiki/index.php?title=conmat.m)
