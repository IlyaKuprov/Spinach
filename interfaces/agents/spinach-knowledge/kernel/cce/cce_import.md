# kernel/cce/cce_import.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_import.m)

## Purpose

Import versioned CCE results, spin tables, tensors and cube samples.

## Behaviour

- MAT/JSON imports recover the model/result envelope; CSV may use the full JSON sidecar or an explicit id,isotope,x,y,z table.
- XYZ coordinates, tensor CSV (isotropic scalar or Cartesian entries), and cube voxel data require explicit units. Canonical imported coordinates are Angstrom and tensor couplings are Hz.
- Optional frame maps input axes to laboratory axes as r_lab=R*r_input and T_lab=R*T_input*R-transpose. Negative-atom-count/orbital cubes are refused.

## Inputs and outputs

```matlab
imported=cce_import(filename,format,units)
```

**Inputs**

- `filename`, `format`, `units` — input path; mat/json/csv/xyz/tensor/cube; external length unit angstrom/nm/m/bohr, tensor unit Hz/kHz/MHz/rad/s/T/mT/G. Field units also require gamma_rad_per_s_T; optional frame is a proper rotation matrix.

**Outputs**

- `imported` — recovered model/result or canonical spin/tensor/cube data with explicit metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_import.m)
- [cce_export.m](cce_export.md)
- [cce_bath.m](cce_bath.md)

Tagged JSON reconstructs empty cells and field-name lists before cell operations, including empty structures, integer arrays, multidimensional zero-sized arrays, and character shapes. MAT, JSON, and CSV-sidecar exports share exact type/shape round-trip tests. Unit fields must belong to the requested physical format; frame and cube density labels remain metadata.
