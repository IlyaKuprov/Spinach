# interfaces/jsonlab-1.5/saveubjson.m

- Signature: `json=saveubjson(rootname,obj,varargin)`

## Purpose

`saveubjson` serializes a MATLAB array, cell array, structure, or class instance as a Universal Binary JSON (UBJSON) binary string.

## Inputs and options

- `rootname`: root-object name. An empty value omits the root name unless `ForceRootName` is enabled; in that case the MATLAB variable name is used.
- `obj`: MATLAB value to serialize.
- Optional filename, options structure, or case-sensitive parameter/value options. The documented options include `FileName`, `ArrayToStruct`, `ParseLogical`, `SingletArray`, `SingletCell`, `ForceRootName`, `JSONP`, and `UnpackHex`.

## Output

`json` is the serialized UBJSON binary string; `FileName` can also direct output to a file.

## Source notes

JSONLab utility by Qianqian Fang; created 2013/08/17. Distributed under the BSD License (see `LICENSE_BSD.txt`).
