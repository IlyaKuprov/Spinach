# interfaces/jsonlab-1.5/struct2jdata.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/struct2jdata.m`
- Signature: `newdata=struct2jdata(data,varargin)`
- Total lines: 96

## Purpose

No descriptive header was found. The best immediate identifier is `newdata=struct2jdata(data,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content

## Implementation structure

- newdata=struct2jdata(data,opt,...)
- convert a JData object (in the form of a struct array) into an array
- authors:Qianqian Fang (q.fang <at> neu.edu)
- input:
- data: a struct array. If data contains JData keywords in the first
- level children, these fields are parsed and regrouped into a
- data object (arrays, trees, graphs etc) based on JData
- specification. The JData keywords are
- "_ArrayType_", "_ArraySize_", "_ArrayData_"
- "_ArrayIsSparse_", "_ArrayIsComplex_"
- opt: (optional) a list of 'Param',value pairs for additional options
- The supported options include

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fieldnames()`, `jsonopt()`, `isstruct()`, `getfield()`, `data()`, `newdata()`, `setfield()`, `jstruct2array()`, `strmatch()`, `cast()`, `double()`, `any()`, `ndata()`, `complex()`, `dim()`, `elseif()`.
