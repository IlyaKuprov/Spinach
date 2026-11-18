# kernel/cce/cce_import.m

## Native role and invocation

Explicit versioned result, spin-table, and Cartesian-tensor import utility. Imported tables are not automatically converted into a native experiment request.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `imported=cce_import(filename,format,units)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_import.m)

## Purpose

Import versioned CCE results, spin tables, tensors, and cube samples.

## Behaviour

- MAT/JSON imports recover the model/result envelope; CSV may use the full JSON sidecar or an explicit id,isotope,x,y,z table.
- XYZ coordinates, tensor CSV (isotropic scalar or Cartesian entries), and cube voxel data require explicit units. Canonical imported coordinates are Angstrom and tensor couplings are Hz.
- The required frame string labels the input frame; optional units.rotation maps input axes to laboratory axes as r_lab=R*r_input and T_lab=R*T_input*R-transpose. Negative-atom-count/orbital cubes are refused.

## Inputs and outputs

**Inputs**

- `filename`, `format`, `units` — input path; mat/json/csv/xyz/tensor/cube; external length unit angstrom/nm/m/bohr, tensor unit Hz/kHz/MHz/rad/s/T/mT/G. Field units also require units.gyro_rad_s_T. units.frame is a required descriptive string for external physical data; units.rotation, if supplied, is a proper rotation matrix.

**Outputs**

- `imported` — recovered model/result or canonical spin/tensor/cube data with explicit metadata.

## References

- [cce_export.m](cce_export.md)

Tagged JSON reconstructs empty cells and field-name lists before cell operations, including empty structures, integer arrays, multidimensional zero-sized arrays, and character shapes. MAT, JSON, and CSV-sidecar exports share exact type/shape round-trip tests. Unit fields must belong to the requested physical format; frame and cube density labels remain metadata.

The retained service uses the native compilation/input boundary: physical models come from cce_compile, and generated baths enter canonical sys/inter through spin_bath.

Native experiment contracts are reconstructed through cce_capabilities; a different capability version is rejected rather than silently substituted. Registered sequence handles and @clexpn deferred descriptors are restored as callable handles. Tagged sparse arrays preserve numeric/logical type and dimensions, character matrices retain column-major layout, and string missingness is retained.

External frame and unit declarations are validated before file parsing; self-describing MAT/JSON envelopes and CSV sidecars remain exempt from external physical-unit declarations. Field-valued tensors retain the explicit gyro_rad_s_T conversion to Hz.
