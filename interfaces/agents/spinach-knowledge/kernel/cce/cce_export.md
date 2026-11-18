# kernel/cce/cce_export.m

## Native role and invocation

Explicit result/model interchange utility, not an automatic clexpn output stage.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `cce_export(filename,model,result)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_export.m)

## Purpose

Export a complete complex CCE result and its model.

## Behaviour

- MAT/JSON store a versioned envelope; JSON preserves numeric type, shape, complex components, and nonfinite values.
- CSV writes top-level numeric result fields as field,index,real,imag rows with column-major linear indices.
- A CSV JSON sidecar at filename.json retains the full result and diagnostics and is the round-trip authority; CSV alone is a human-readable view.

## Inputs and outputs

**Inputs**

- `filename` — .mat, .json, or .csv destination.
- `model`, `result` — scalar structures to serialise.

**Outputs**

- No return value; writes the requested file (and JSON sidecar for CSV).

## References

- [cce_import.m](cce_import.md)

Character tags retain MATLAB array shape, including nonstandard empty character dimensions. The importer remains compatible with earlier character tags without a shape field.

JSON and CSV sidecars preserve complete native diagnostics: registered capability contracts are stored by experiment name and capability version, then reconstructed through cce_capabilities on import. Registered sequence handles and @clexpn deferred descriptors are portable; arbitrary handles are rejected. Sparse arrays retain storage type, shape, and complex values; character matrices and missing string elements retain their original values.
