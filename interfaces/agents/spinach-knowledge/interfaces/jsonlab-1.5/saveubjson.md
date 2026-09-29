# interfaces/jsonlab-1.5/saveubjson.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/jsonlab-1.5/saveubjson.m) · [JSONLab project page](http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)

**Call:** `json=saveubjson(rootname,obj,varargin)`. The implementation also accepts `saveubjson(obj)`; the one-input form derives a root name from the input variable, or uses `root` when that value is unnamed. With a root name, use `saveubjson(rootname,obj,filename)`, an options structure, or name/value pairs. A sole character optional argument is interpreted as a filename; `FileName` is the corresponding option.

The accepted MATLAB data families are numeric/logical arrays, character data, cells, structs/struct arrays, and class instances; class instances are serialised from their visible properties. An empty root name omits the root wrapper unless `ForceRootName` is enabled. The returned `json` is a character vector containing UBJSON bytes; `FileName` writes it in binary mode. Numeric arrays above two dimensions, sparse arrays, complex arrays, and arrays selected by `ArrayToStruct` use JData metadata; higher-dimensional cells and struct arrays are reshaped to two dimensions. The serializer does not assign or convert physical units.

Options read by this implementation include `FileName`, `ArrayToStruct` (default 0), `SingletArray` (default 0), `SingletCell` (default 1), `ForceRootName` (default 0), `JSONP`, `Inf`, `NaN`, and `UnpackHex` (default 1). `ArrayToStruct` stores array type and size with data; sparse arrays use index/value triplets and complex arrays carry real/imaginary components and a complex marker. `JSONP`, when nonempty, wraps the serialised result in the named call form. The source code does not consult several options mentioned in copied header text, including `FloatFormat`, `ArrayIndent`, and `ParseLogical`; they are not described here as supported controls.

The source example builds a mesh struct with coordinate/connectivity arrays and `SpecialData=[NaN,Inf,-Inf]`, then calls `saveubjson('jsonmesh',jsonmesh)` or writes `meshdata.ubj`. The example is source documentation, not a run result. Non-finite-value encodings and downstream interpretation should be checked against the consumer's UBJSON implementation.

The vendored source identifies Qianqian Fang as author and gives a BSD license (see `LICENSE_BSD.txt`).
