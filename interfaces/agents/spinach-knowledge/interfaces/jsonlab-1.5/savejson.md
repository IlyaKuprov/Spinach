# interfaces/jsonlab-1.5/savejson.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/savejson.m`
- Signature: `json=savejson(rootname,obj,varargin)`
- Total lines: 552

## Purpose

No descriptive header was found. The best immediate identifier is `json=savejson(rootname,obj,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file also defines local helper function(s): `obj2json()`, `cell2json()`, `struct2json()`, `str2json()`, `mat2json()`, `matlabobject2json()`, `matdata2json()`, `checkname()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 2-92: json=savejson(rootname,obj,filename) json=savejson(rootname,obj,opt) json=savejson(rootname,obj,'param1',value1,'param2',value2,...) convert a MATLAB object (cell, struct or array) into a JSON (JavaScript Object Notation) string created on 2011/09/09 $Id$ input: rootname: the name of the root-object, when set to '', the root name is ignored, however, when opt.ForceRootName is set to 1 (see below), the MATLAB variable name will be used as the root name. obj: a MATLAB object (array, cell, cell array, struct, struct array, class instance). filename: a string for the file name to save the output JSON data. opt: a struct for additional options, ignore to use default values. opt can have the following fields (first in [.|.] is the default) opt.FileName [''|string]: a file name to save the output JSON data opt.FloatFormat ['%.10g'|string]: format to show each numeric element of a 1D/2D array; opt.ArrayIndent [1|0]: if 1, output explicit data array with precedent indentation; if 0, no indentation opt.ArrayToStruct[0|1]: when set to 0, savejson outputs 1D/2D array in JSON array format; if sets to 1, an array will be shown as a struct with fields "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for sparse arrays, the non-zero elements will be saved to _ArrayData_ field in triplet-format i.e. (ix,iy,val) and "_ArrayIsSparse_" will be added with a value of 1; for a complex array, the _ArrayData_ array will include two columns (4 for sparse) to record the real and imaginary parts, and also "_ArrayIsComplex_":1 is added. opt.ParseLogical [0|1]: if this is set to 1, logical array elem will use true/false rather than 1/0. opt.SingletArray [0|1]: if this is set to 1, arrays with a single numerical element will be shown without a square bracket, unless it is the root object; if 0, square brackets are forced for any numerical arrays. opt.SingletCell [1|0]: if 1, always enclose a cell with "[]" even it has only one element; if 0, brackets are ignored when a cell has only 1 element. opt.ForceRootName [0|1]: when set to 1 and rootname is empty, savejson will use the name of the passed obj variable as the root object name; if obj is an expression and does not have a name, 'root' will be used; if this is set to 0 and rootname is empty, the root level will be merged down to the lower level. opt.Inf ['"$1_Inf_"'|string]: a customized regular expression pattern to represent +/-Inf. The matched pattern is '([-+]*)Inf' and $1 represents the sign. For those who want to use 1e999 to represent Inf, they can set opt.Inf to '$11e999' opt.NaN ['"_NaN_"'|string]: a customized regular expression pattern to represent NaN opt.JSONP [''|string]: to generate a JSONP output (JSON with padding), for example, if opt.JSONP='foo', the JSON data is wrapped inside a function call as 'foo(...);' opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson back to the string form opt.SaveBinary [0|1]: 1 -save the JSON file in binary mode; 0 -text mode. opt.Compact [0|1]: 1-out compact JSON format (remove all newlines and tabs) opt can be replaced by a list of ('param',value) pairs. The param string is equivallent to a field in opt and is case sensitive. output: examples: jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 'MeshTetra',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],... 'MeshTri',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;... 2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],... 'MeshCreator','FangQ','MeshTitle','T6 Cube',... 'SpecialData',[nan, inf, -inf]); savejson('jmesh',jsonmesh) savejson('',jsonmesh,'ArrayIndent',0,'FloatFormat','\t%.5g') license: BSD License, see LICENSE_BSD.txt files for details; implemented by `if(nargin==1)`.
- Lines 152-153: save to a file if FileName is set, suggested by Patrick Rapin; implemented by `filename=jsonopt('FileName','',opt)`.

### Key state/data transformations

- Lines 93: computes `varname` using `varname=inputname(1)`.
- Lines 94: computes `obj` using `obj=rootname`.
- Lines 98: computes `rootname` using `rootname=varname`.
- Lines 103: computes `opt` using `opt=struct('filename',varargin{1})`.
- Lines 107: computes `opt.IsOctave` using `opt.IsOctave=exist('OCTAVE_VERSION','builtin')`.
- Lines 111: computes `opt.singletarray` using `opt.singletarray=not(opt.norowbracket)`.
- Lines 114: computes `rootisarray` using `rootisarray=0`.
- Lines 115: computes `rootlevel` using `rootlevel=1`.
- Lines 116: computes `forceroot` using `forceroot=jsonopt('ForceRootName',0,opt)`.
- Lines 130: computes `whitespaces` using `whitespaces=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'))`.
- Lines 135: computes `opt.whitespaces_` using `opt.whitespaces_=whitespaces`.
- Lines 138: computes `nl` using `nl=whitespaces.newline`.
- Lines 140: computes `json` using `json=obj2json(rootname,obj,rootlevel,opt)`.
- Lines 147: computes `jsonp` using `jsonp=jsonopt('JSONP','',opt)`.
- Lines 153: computes `filename` using `filename=jsonopt('FileName','',opt)`.
- Lines 156: computes `fid` using `fid = fopen(filename, 'wb')`.

### Local helper functions

- Line 166: `obj2json()` — `function txt=obj2json(name,item,level,varargin)`.
  - Representative operation: `if(iscell(item))`.
  - Representative operation: `txt=cell2json(name,item,level,varargin{:})`.
- Line 181: `cell2json()` — `function txt=cell2json(name,item,level,varargin)`.
  - Representative operation: `txt={}`.
  - Representative operation: `if(~iscell(item))`.
- Line 235: `struct2json()` — `function txt=struct2json(name,item,level,varargin)`.
  - Representative operation: `txt={}`.
  - Representative operation: `if(~isstruct(item))`.
- Line 311: `str2json()` — `function txt=str2json(name,item,level,varargin)`.
  - Representative operation: `txt={}`.
  - Representative operation: `if(~ischar(item))`.
- Line 356: `mat2json()` — `function txt=mat2json(name,item,level,varargin)`.
  - Representative operation: `if(~isnumeric(item) && ~islogical(item))`.
  - Representative operation: `error('input is not an array')`.
- Line 434: `matlabobject2json()` — `function txt=matlabobject2json(name,item,level,varargin)`. "st = struct(item);" would produce an inmutable warning, because it make the protected and private properties visible. Instead we get the
  - Representative operation: `if numel(item) == 0`.
  - Representative operation: `st = struct()`.
- Line 451: `matdata2json()` — `function txt=matdata2json(mat,level,varargin)`.
  - Representative operation: `ws=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'))`.
  - Representative operation: `ws=jsonopt('whitespaces_',ws,varargin{:})`.
- Line 502: `checkname()` — `function newname=checkname(name,varargin)`.
  - Representative operation: `isunpack=jsonopt('UnpackHex',1,varargin{:})`.
  - Representative operation: `newname=name`.

## Implementation structure

- json=savejson(rootname,obj,filename)
- json=savejson(rootname,obj,opt)
- json=savejson(rootname,obj,'param1',value1,'param2',value2,...)
- convert a MATLAB object (cell, struct or array) into a JSON (JavaScript
- Object Notation) string
- author: Qianqian Fang (q.fang <at> neu.edu)
- created on 2011/09/09
- $Id$
- input:
- rootname: the name of the root-object, when set to '', the root name
- is ignored, however, when opt.ForceRootName is set to 1 (see below),
- the MATLAB variable name will be used as the root name.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `inputname()`, `ischar()`, `varargin2struct()`, `exist()`, `isfield()`, `not()`, `jsonopt()`, `islogical()`, `isstruct()`, `iscell()`, `isobject()`, `obj2json()`, `fopen()`, `fwrite()`, `fclose()`, `cell2json()`.
