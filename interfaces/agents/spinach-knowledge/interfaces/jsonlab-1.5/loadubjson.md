# interfaces/jsonlab-1.5/loadubjson.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/loadubjson.m`
- Signature: `data = loadubjson(fname,varargin)`
- Total lines: 457

## Purpose

No descriptive header was found. The best immediate identifier is `data = loadubjson(fname,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content

- The file also defines local helper function(s): `parse_object()`, `elem_info()`, `parse_block()`, `parse_array()`, `parse_char()`, `parse_name()`, `parseStr()`, `parse_number()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 2-52: data=loadubjson(fname,opt) data=loadubjson(fname,'param1',value1,'param2',value2,...) parse a JSON (JavaScript Object Notation) file or string created on 2013/08/01 $Id$ input: fname: input file name, if fname contains "{}" or "[]", fname will be interpreted as a UBJSON string opt: a struct to store parsing options, opt can be replaced by a list of ('param',value) pairs -the param string is equivallent to a field in opt. opt can have the following fields (first in [.|.] is the default) opt.SimplifyCell [0|1]: if set to 1, loadubjson will call cell2mat for each element of the JSON data, and group arrays based on the cell2mat rules. opt.IntEndian [B|L]: specify the endianness of the integer fields in the UBJSON input data. B -Big-Endian format for integers (as required in the UBJSON specification); L -input integer fields are in Little-Endian order. opt.NameIsString [0|1]: for UBJSON Specification Draft 8 or earlier versions (JSONLab 1.0 final or earlier), the "name" tag is treated as a string. To load these UBJSON data, you need to manually set this flag to 1. output: dat: a cell array, where {...} blocks are converted into cell arrays, and [...] are converted to arrays examples: obj=struct('string','value','array',[1 2 3]); ubjdata=saveubjson('obj',obj); dat=loadubjson(ubjdata) dat=loadubjson(['examples' filesep 'example1.ubj']) dat=loadubjson(['examples' filesep 'example1.ubj'],'SimplifyCell',1) license: BSD License, see LICENSE_BSD.txt files for details; implemented by `global pos inStr len esc index_esc len_esc isoct arraytoken fileendian systemendian`.
- Lines 71-72: String delimiters and escape chars identified to improve speed:; implemented by `esc = find(inStr=='"' | inStr=='\' )`.

### Control flow inferred from the code

- Line 80: `while` loop over `pos <= len`.

### Key state/data transformations

- Lines 55: computes `string` using `string=fname`.
- Lines 57: computes `fid` using `fid = fopen(fname,'rb')`.
- Lines 64: computes `pos` using `pos = 1; len = length(string); inStr = string`.
- Lines 65: computes `isoct` using `isoct=exist('OCTAVE_VERSION','builtin')`.
- Lines 67: computes `jstr` using `jstr=regexprep(inStr,'\\\\',' ')`.
- Lines 68: computes `escquote` using `escquote=regexp(jstr,'\\"')`.
- Lines 69: computes `arraytoken` using `arraytoken=sort([arraytoken escquote])`.
- Lines 73: computes `index_esc` using `index_esc = 1; len_esc = length(esc)`.
- Lines 75: computes `opt` using `opt=varargin2struct(varargin{:})`.
- Lines 76: computes `fileendian` using `fileendian=upper(jsonopt('IntEndian','B',opt))`.
- Lines 77: computes `[os,maxelem,systemendian]` using `[os,maxelem,systemendian]=computer`.
- Lines 79: computes `jsoncount` using `jsoncount=1`.
- Lines 83: computes `data{jsoncount}` using `data{jsoncount} = parse_object(opt)`.
- Lines 94: computes `data` using `data=data{1}`.

### Local helper functions

- Line 98: `parse_object()` — `function object = parse_object(varargin)`.
  - Representative operation: `parse_char('{')`.
  - Representative operation: `object = []`.
- Line 140: `elem_info()` — `function [cid,len]=elem_info(type)`. %-------------------------------------------------------------------------
  - Representative operation: `id=strfind('iUIlLdD',type)`.
  - Representative operation: `dataclass={'int8','uint8','int16','int32','int64','single','double'}`.
- Line 153: `parse_block()` — `function [data, adv]=parse_block(type,count,varargin)`.
  - Representative operation: `global pos inStr isoct fileendian systemendian`.
  - Representative operation: `[cid,len]=elem_info(type)`.
- Line 172: `parse_array()` — `function object = parse_array(varargin)`.
  - Representative operation: `global pos inStr`.
  - Representative operation: `parse_char('[')`.
- Line 238: `parse_char()` — `function parse_char(c)`. %-------------------------------------------------------------------------
  - Representative operation: `global pos inStr len`.
  - Representative operation: `skip_whitespace`.
- Line 250: `function c = next_char()` — `function c = next_char`. %-------------------------------------------------------------------------
  - Representative operation: `global pos inStr len`.
  - Representative operation: `skip_whitespace`.
- Line 261: `function skip_whitespace()` — `function skip_whitespace`. %-------------------------------------------------------------------------
  - Representative operation: `global pos inStr len`.
  - Representative operation: `while pos <= len && isspace(inStr(pos))`.
- Line 268: `parse_name()` — `function str = parse_name(varargin)`. %-------------------------------------------------------------------------
  - Representative operation: `global pos inStr`.
  - Representative operation: `bytelen=double(parse_number())`.

## Implementation structure

- data=loadubjson(fname,opt)
- data=loadubjson(fname,'param1',value1,'param2',value2,...)
- parse a JSON (JavaScript Object Notation) file or string
- authors:Qianqian Fang (q.fang <at> neu.edu)
- created on 2013/08/01
- $Id$
- input:
- fname: input file name, if fname contains "{}" or "[]", fname
- will be interpreted as a UBJSON string
- opt: a struct to store parsing options, opt can be replaced by
- a list of ('param',value) pairs -the param string is equivallent
- to a field in opt. opt can have the following

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `regexp()`, `elseif()`, `exist()`, `fopen()`, `fread()`, `fclose()`, `regexprep()`, `varargin2struct()`, `upper()`, `jsonopt()`, `parse_object()`, `parse_array()`, `error_pos()`, `iscell()`, `parse_char()`, `inStr()`.
