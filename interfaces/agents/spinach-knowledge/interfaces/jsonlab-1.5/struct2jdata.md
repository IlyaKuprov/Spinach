# interfaces/jsonlab-1.5/struct2jdata.m

- Signature: `newdata=struct2jdata(data,varargin)`

## Purpose

Converts a JData object represented by a struct array into an array, regrouping its first-level JData keyword fields according to the JData specification.

## Input and options

- `data`: struct array. Recognized array fields include `_ArrayType_`, `_ArraySize_`, `_ArrayData_`, `_ArrayIsSparse_`, and `_ArrayIsComplex_`.
- `Recursive` (optional parameter): 1 applies conversion to every child; 0 disables recursive conversion. Other optional parameter/value pairs are accepted.

## Output

`newdata` contains the converted data when the input has a JData structure; otherwise the input is returned unchanged.

## Source notes

JSONLab toolbox utility by Qianqian Fang (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab); the source specifies a BSD license.
