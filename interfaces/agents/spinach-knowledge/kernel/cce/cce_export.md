# kernel/cce/cce_export.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_export.m)

## Purpose

Export a complete complex CCE result and its model.

## Behaviour

- MAT/JSON store a versioned envelope; JSON preserves numeric type, shape, complex components and nonfinite values.
- CSV writes numeric vector/matrix leaves of result as leaf,index1,index2,real,imag rows.
- A CSV JSON sidecar at filename.json retains the full result and diagnostics and is the round-trip authority; CSV alone is a human-readable view.

## Inputs and outputs

```matlab
cce_export(filename,model,result)
```

**Inputs**

- `filename` — .mat, .json or .csv destination.
- `model`, `result` — scalar structures to serialize.

**Outputs**

- No return value; writes the requested file (and JSON sidecar for CSV).

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_export.m)
- [cce_import.m](cce_import.md)
