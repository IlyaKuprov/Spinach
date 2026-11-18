# interfaces/jsonlab-1.5/struct2jdata.m

This source is absent from the current distribution. The following historical description is retained for reference, not as an available runnable interface. No replacement is inferred from the old signature.

`interfaces/jsonlab-1.5/struct2jdata.m` (source removed) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `newdata=struct2jdata(data,varargin)` converts a JSONLab JData struct representation back to MATLAB data. `data` is expected to be a struct or struct array. The implementation checks the MATLAB-side field names `x0x5F_ArrayType_` and `x0x5F_ArrayData_`; optional metadata fields are `x0x5F_ArraySize_`, `x0x5F_ArrayIsSparse_`, and `x0x5F_ArrayIsComplex_`.

For an array record it casts `_ArrayData_` to the declared array type. For dense records without an `_ArrayIsSparse_` field, size metadata triggers `reshape`, and two-column complex data are combined into complex values first. Sparse records with a true sparse flag are rebuilt from index/value data, using declared dimensions when supplied or inferring them otherwise; complex sparse records combine real and imaginary value columns when the encoded width matches the expected representation. A present but false sparse flag bypasses the dense reshape branch; the function therefore depends on the metadata representation, not just the flag value. If there is no recognised array-type/data pair, the original structure is retained apart from any requested child conversion. `Recursive` defaults to 0; when set to 1, the implementation passes each struct-valued child field to `jstruct2array` before considering the current record. It does not call `struct2jdata` recursively, and it does not forward the options to `jstruct2array`; deeper traversal is therefore not specified by this function alone.

For recognised array records, a one-element input returns the converted value directly; a multiple-element input returns a cell column of converted values and iterates over `length(data)` (not `numel(data)`), so the intended input is a vector-like struct array. Type casts, shape reconstruction, and sparse construction rely on consistent JData metadata; the function does not add domain-specific validation or unit conversion. Recognised array records are replaced by the reconstructed array, so unrelated payload fields are not carried into that array.

The source's example describes a sparse double array with `_ArraySize_=[2 3]`, `_ArrayIsSparse_=1`, and `_ArrayData_` as the encoded data field. The source example is documentation, not a claim that it was executed. The source uses JSONLab's escaped MATLAB field names in its lookups; provide data in the representation produced by the corresponding JSONLab loader.

The vendored source identifies Qianqian Fang as author and gives a BSD licence (see `LICENSE_BSD.txt`).
