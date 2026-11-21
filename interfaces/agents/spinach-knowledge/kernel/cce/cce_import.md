# kernel/cce/cce_import.m

## Native role and invocation

Explicit import of spin tables, Cartesian tensors, and cube density samples. Imported tables are not automatically converted into a native experiment request.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `imported=cce_import(filename,format,units)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_import.m)

## Purpose

Import external spin tables, tensors, and cube samples into canonical physical units.

## Behaviour

- Spin CSV requires the columns `id,isotope,x,y,z`; it is read directly and retains physical IDs and isotope labels.
- XYZ coordinates, tensor CSV (isotropic scalar or Cartesian entries), and cube voxel data require explicit units. Canonical imported coordinates are Angstrom and tensor couplings are Hz.
- The required frame string labels the input frame; optional units.rotation maps input axes to laboratory axes as r_lab=R*r_input and T_lab=R*T_input*R-transpose. Negative-atom-count/orbital cubes are refused.

## Inputs and outputs

**Inputs**

- `filename`, `format`, `units` — input path; csv/xyz/tensor/cube; external length unit angstrom/nm/m/bohr, tensor unit Hz/kHz/MHz/rad/s/T/mT/G. Field units also require units.gyro_rad_s_T. units.frame is a required descriptive string for external physical data; units.rotation, if supplied, is a proper rotation matrix.

**Outputs**

- `imported` — canonical spin/tensor/cube data with explicit metadata.

Unit fields must belong to the requested physical format. Frame and cube density labels remain metadata. Frame and unit declarations are validated before file parsing; field-valued tensors use the explicit `gyro_rad_s_T` conversion to Hz.

Physical models come from `cce_compile`, and generated baths enter canonical `sys`/`inter` through `spin_bath`.
