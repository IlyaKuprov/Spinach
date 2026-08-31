# interfaces/jsonlab-1.5/loadjson.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/loadjson.m`
- Signature: `data = loadjson(fname,varargin)`
- Total lines: 503

## Purpose

No descriptive header was found. The best immediate identifier is `data = loadjson(fname,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content

- The file also defines local helper function(s): `parse_object()`, `parse_array()`, `skip_whitespace()`, `next_char()`, `parseStr()`, `parse_number()`, `parse_value()`, `max()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 2-64: data=loadjson(fname,opt) data=loadjson(fname,'param1',value1,'param2',value2,...) parse a JSON (JavaScript Object Notation) file or string created on 2011/09/09, including previous works from created on 2009/11/02 created on 2009/03/22 Joel Feenstra: created on 2008/07/03 $Id$ input: fname: input file name, if fname contains "{}" or "[]", fname will be interpreted as a JSON string opt: a struct to store parsing options, opt can be replaced by a list of ('param',value) pairs -the param string is equivallent to a field in opt. opt can have the following fields (first in [.|.] is the default) opt.SimplifyCell [0|1]: if set to 1, loadjson will call cell2mat for each element of the JSON data, and group arrays based on the cell2mat rules. opt.FastArrayParser [1|0 or integer]: if set to 1, use a speed-optimized array parser when loading an array object. The fast array parser may collapse block arrays into a single large array similar to rules defined in cell2mat; 0 to use a legacy parser; if set to a larger-than-1 value, this option will specify the minimum dimension to enable the fast array parser. For example, if the input is a 3D array, setting FastArrayParser to 1 will return a 3D array; setting to 2 will return a cell array of 2D arrays; setting to 3 will return to a 2D cell array of 1D vectors; setting to 4 will return a 3D cell array. opt.ShowProgress [0|1]: if set to 1, loadjson displays a progress bar. output: dat: a cell array, where {...} blocks are converted into cell arrays, and [...] are converted to arrays examples: dat=loadjson('{"obj":{"string":"value","array":[1,2,3]}}') dat=loadjson(['examples' filesep 'example1.json']) dat=loadjson(['examples' filesep 'example1.json'],'SimplifyCell',1) license: BSD License, see LICENSE_BSD.txt files for details; implemented by `global pos index_esc isoct arraytoken`.
- Lines 89-90: String delimiters and escape chars identified to improve speed:; implemented by `esc = find(inStr=='"' | inStr=='\' )`.

### Control flow inferred from the code

- Line 99: `while` loop over `pos <= len`.

### Key state/data transformations

- Lines 67: computes `string` using `string=fname`.
- Lines 82: computes `pos` using `pos = 1; len = length(string); inStr = string`.
- Lines 83: computes `isoct` using `isoct=exist('OCTAVE_VERSION','builtin')`.
- Lines 85: computes `jstr` using `jstr=regexprep(inStr,'\\\\',' ')`.
- Lines 86: computes `escquote` using `escquote=regexp(jstr,'\\"')`.
- Lines 87: computes `arraytoken` using `arraytoken=sort([arraytoken escquote])`.
- Lines 91: computes `index_esc` using `index_esc = 1`.
- Lines 93: computes `opt` using `opt=varargin2struct(varargin{:})`.
- Lines 96: computes `opt.progressbar_` using `opt.progressbar_=waitbar(0,'loading ')`.
- Lines 98: computes `jsoncount` using `jsoncount=1`.
- Lines 102: computes `data{jsoncount}` using `data{jsoncount} = parse_object(inStr, esc, opt)`.
- Lines 113: computes `data` using `data=data{1}`.

### Local helper functions

- Line 121: `parse_object()` — `function object = parse_object(inStr, esc, varargin)`.
  - Representative operation: `parse_char(inStr, '{')`.
  - Representative operation: `object = []`.
- Line 146: `parse_array()` — `function object = parse_array(inStr, esc, varargin)`.
  - Representative operation: `global pos isoct`.
  - Representative operation: `parse_char(inStr, '[')`.
- Line 251: `parse_char()` — `function parse_char(inStr, c)`. %-------------------------------------------------------------------------
  - Representative operation: `global pos`.
  - Representative operation: `pos=skip_whitespace(pos, inStr)`.
- Line 263: `next_char()` — `function c = next_char(inStr)`. %-------------------------------------------------------------------------
  - Representative operation: `global pos`.
  - Representative operation: `pos=skip_whitespace(pos, inStr)`.
- Line 274: `skip_whitespace()` — `function newpos=skip_whitespace(pos, inStr)`. %-------------------------------------------------------------------------
  - Representative operation: `newpos=pos`.
  - Representative operation: `while newpos <= length(inStr) && isspace(inStr(newpos))`.
- Line 281: `parseStr()` — `function str = parseStr(inStr, esc, varargin)`. len, ns = length(inStr), keyboard
  - Representative operation: `global pos index_esc`.
  - Representative operation: `if inStr(pos) ~= '"'`.
- Line 345: `parse_number()` — `function num = parse_number(inStr, varargin)`.
  - Representative operation: `global pos isoct`.
  - Representative operation: `currstr=inStr(pos:min(pos+30,end))`.
- Line 362: `parse_value()` — `function val = parse_value(inStr, esc, varargin)`.
  - Representative operation: `global pos`.
  - Representative operation: `len=length(inStr)`.

## Implementation structure

- data=loadjson(fname,opt)
- data=loadjson(fname,'param1',value1,'param2',value2,...)
- parse a JSON (JavaScript Object Notation) file or string
- authors:Qianqian Fang (q.fang <at> neu.edu)
- created on 2011/09/09, including previous works from
- Nedialko Krouchev: http://www.mathworks.com/matlabcentral/fileexchange/25713
- created on 2009/11/02
- Francois Glineur: http://www.mathworks.com/matlabcentral/fileexchange/23393
- created on 2009/03/22
- Joel Feenstra:
- created on 2008/07/03
- $Id$

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `regexp()`, `elseif()`, `exist()`, `fileread()`, `urlread()`, `fullfile()`, `regexprep()`, `varargin2struct()`, `jsonopt()`, `waitbar()`, `next_char()`, `parse_object()`, `parse_array()`, `error_pos()`, `iscell()`, `isfield()`.
