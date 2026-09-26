# interfaces/jsonlab-1.5/varargin2struct.m

- Signature: `opt=varargin2struct(varargin)`

## Purpose

No descriptive header was found. The best immediate identifier is `opt=varargin2struct(varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content

## Implementation structure

- opt=varargin2struct('param1',value1,'param2',value2,...)
- opt=varargin2struct(...,optstruct,...)
- convert a series of input parameters into a structure
- authors:Qianqian Fang (q.fang <at> neu.edu)
- date: 2012/12/22
- input:
- 'param', value: the input parameters should be pairs of a string and a value
- optstruct: if a parameter is a struct, the fields will be merged to the output struct
- output:
- opt: a struct where opt.param1=value1, opt.param2=value2 ...
- license:
- BSD License, see LICENSE_BSD.txt files for details
