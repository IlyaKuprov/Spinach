# interfaces/jsonlab-1.5/struct2jdata.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/struct2jdata.m) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `newdata=struct2jdata(data,varargin)` converts a JSONLab JData struct representation back to MATLAB data. `data` is expected to be a struct or struct array. The implementation checks the MATLAB-side field names `x0x5F_ArrayType_` and `x0x5F_ArrayData_`; optional metadata fields are `x0x5F_ArraySize_`, `x0x5F_ArrayIsSparse_`, and `x0x5F_ArrayIsComplex_`.

For an array record it casts `_ArrayData_` to the declared array type. With size metadata it reconstructs the shape using `reshape`; complex records combine real/imaginary columns. Sparse records are rebuilt as sparse matrices from index/value data and the declared dimensions. If there is no recognised array-type/data pair, the current value is retained. `Recursive` defaults to 0; when set to 1, struct-valued child fields are traversed depth-first before the current record is converted.

For a one-element input, the converted value is returned directly. For inputs with multiple elements, the implementation returns a cell column of converted values and iterates over `length(data)` (not `numel(data)`), so the intended input is a vector-like struct array. Type casts, shape reconstruction, and sparse construction rely on consistent JData metadata; the function does not add domain-specific validation or unit conversion. Units, if present in payload fields, remain payload data.

The source's example describes a sparse double array with `_ArraySize_=[2 3]`, `_ArrayIsSparse_=1`, and `_ArrayData_` as the encoded data field. The source example is documentation, not a claim that it was executed. The source uses JSONLab's escaped MATLAB field names in its lookups; provide data in the representation produced by the corresponding JSONLab loader.

The vendored source identifies Qianqian Fang as author and gives a BSD license (see `LICENSE_BSD.txt`).
