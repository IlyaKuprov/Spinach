# interfaces/jsonlab-1.5/savejson.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/savejson.m) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `json=savejson(rootname,obj,varargin)`. The implementation also accepts `savejson(obj)`; the one-input form uses the input variable name as the root name, falling back to `root` for an unnamed expression. With a root name, use `savejson(rootname,obj,filename)`, an options structure, or name/value pairs. A sole character optional argument is treated as the filename; `FileName` is the equivalent option.

The serializer accepts MATLAB numeric and logical arrays, character data, cells, structs/struct arrays, and class instances; class instances are serialised from their visible properties. An empty `rootname` omits the root wrapper by default; `ForceRootName` requests a root (using the object variable name when available). The result is a JSON character vector. Numeric arrays above two dimensions, sparse arrays, complex arrays, and arrays selected by `ArrayToStruct` use JData metadata; higher-dimensional cells and struct arrays are reshaped to a two-dimensional layout. Singleton brackets depend on `SingletArray` and `SingletCell`. No physical units are inferred, converted, or attached.

Options documented and used by this implementation include:

- `FloatFormat` (default `%.10g`) and `ArrayIndent` (default 1) control numeric formatting and array indentation; `Compact` removes formatting whitespace.
- `ArrayToStruct` (default 0) encodes arrays with JData fields such as `_ArrayType_`, `_ArraySize_`, and `_ArrayData_`. Sparse data use index/value triplets plus `_ArrayIsSparse_`; complex data include real/imaginary components and `_ArrayIsComplex_`.
- `ParseLogical` (default 0) selects `true`/`false` instead of `1`/`0`. `SingletArray` (default 0) controls brackets for singleton numeric arrays, and `SingletCell` (default 1) controls brackets for one-element cells.
- `Inf` and `NaN` choose replacement patterns for non-finite numeric values; defaults serialise them as strings such as `"_Inf_"` and `"_NaN_"`. `JSONP` optionally wraps the output as a callback expression.
- `UnpackHex` (default 1) controls conversion of escaped hexadecimal text. `SaveBinary` (default 0) selects binary or text file writing when `FileName` is set.

The source example serialises a mesh struct containing numeric coordinate/connectivity arrays and `SpecialData=[NaN,Inf,-Inf]`; it also demonstrates `savejson('',jsonmesh,'ArrayIndent',0,'FloatFormat','\t%.5g')`. These are serialisation examples only, not claims of a run. The function does not validate application-specific units or meanings, and JSON has no native NaN/Inf values, so use the replacement options if downstream consumers need a different convention.

The vendored source identifies Qianqian Fang as author and gives a BSD license (see `LICENSE_BSD.txt`).
