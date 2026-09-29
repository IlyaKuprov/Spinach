# interfaces/mestrenova/s2json.m

- Signature: `s2json(file_name,sys,inter,parameters,fid)` (all five inputs are required; the function defines no defaults).

## Purpose and data contract

Serialises the supplied Spinach structures and FID through JSONLab's `savejson` for import into MestreNova. The arguments are placed in a MATLAB structure with fields `sys`, `inter`, `parameters`, and `fid`; `savejson('spinach',spinach,file_name)` writes that structure under the JSON root key `spinach`. The function returns no MATLAB output and performs no unit conversion or FID reshaping.

- `file_name`: character array passed to JSONLab as the destination. The guard checks only `ischar`; it does not check that the path is writable or that its directory exists.
- `sys`, `inter`, `parameters`: MATLAB structures. The implementation checks only that each is a structure; it does not validate fields, array shape, or numerical values.
- `fid`: either numeric data of any shape or a structure. The documented Fourier-transform-only case is a complex matrix. For States quadrature data such as NOESY, supply a structure containing `fid.cos` and `fid.sin` matrices. The code accepts any structure here and does not check those fields or matrix dimensions.

The serialised content retains the supplied parameter and signal representation; this function does not assign physical units. Invalid argument types raise the explicit input errors; serialisation and file-write behaviour are delegated to `savejson`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/mestrenova/s2json.m)
- [Spin Dynamics Wiki: s2json.m](https://spindynamics.org/wiki/index.php?title=s2json.m)
