# interfaces/orca/ocparse.m

- Signature: `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)` (both inputs are required; `pad_factor` has no default).

## Purpose and input format

Parses an ORCA cube file in the documented 3-D simple format and returns an absolute, numerically normalised density array and grid metrics. It imports the file with four header lines and interprets the text metadata as grid point counts, corner coordinates, and three grid-step components. The numeric cube data are reshaped and permuted to `[X Y Z]` order.

- `filename`: nonempty character string naming an existing file. The guard checks existence, but does not validate that the file is a readable ORCA cube or that its header and data are internally consistent.
- `pad_factor`: nonnegative real integer scalar. No padding is added at zero. At a positive value, each side of each dimension receives `pad_factor` times that dimension's point count in zeros, so the padded array shape is `(1+2*pad_factor)*[nx ny nz]`.

## Outputs and units

- `density`: padded array in `[X Y Z]` axis order. The parser takes `abs(A.data)` and normalises the unpadded array by three trapezoidal integrations multiplied by `dx*dy*dz`; zero padding is then applied. The code does not separately guard against a zero or non-finite normalisation integral.
- `ext`: six coordinate values in the order `[xmin xmax ymin ymax zmin zmax]`. Before padding, each pair is computed from the corresponding corner coordinate, point count, and grid step. Padding expands the pair by `npoints*step*pad_factor` on each side.
- `dx`, `dy`, `dz`: the three grid-step components read from the cube metadata, in Angstrom. The coordinate extents are also in Angstrom. Density units are not separately declared by the function; its normalisation uses the numerical grid-step product.

The only explicit guards are for filename type/non-emptiness, file existence, and pad-factor type/reality/scalar/non-negativity/integrality. Malformed cube content or incompatible dimensions are left to the importer and subsequent MATLAB operations to reject.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/orca/ocparse.m)
- [Spin Dynamics Wiki: ocparse.m](https://spindynamics.org/wiki/index.php?title=ocparse.m)
- Authors credited by the source: Ilya Kuprov, Elizaveta Suturina, and Petra Pikulova.
