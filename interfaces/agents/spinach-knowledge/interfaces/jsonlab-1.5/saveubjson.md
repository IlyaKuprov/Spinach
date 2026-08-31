# interfaces/jsonlab-1.5/saveubjson.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/saveubjson.m`
- Signature: `json=saveubjson(rootname,obj,varargin)`
- Total lines: 562

## Purpose

No descriptive header was found. The best immediate identifier is `json=saveubjson(rootname,obj,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file also defines local helper function(s): `obj2ubjson()`, `cell2ubjson()`, `struct2ubjson()`, `str2ubjson()`, `mat2ubjson()`, `matlabobject2ubjson()`, `matdata2ubjson()`, `checkname()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 2-80: json=saveubjson(rootname,obj,filename) json=saveubjson(rootname,obj,opt) json=saveubjson(rootname,obj,'param1',value1,'param2',value2,...) convert a MATLAB object (cell, struct or array) into a Universal Binary JSON (UBJSON) binary string created on 2013/08/17 $Id$ input: rootname: the name of the root-object, when set to '', the root name is ignored, however, when opt.ForceRootName is set to 1 (see below), the MATLAB variable name will be used as the root name. obj: a MATLAB object (array, cell, cell array, struct, struct array, class instance) filename: a string for the file name to save the output UBJSON data opt: a struct for additional options, ignore to use default values. opt can have the following fields (first in [.|.] is the default) opt.FileName [''|string]: a file name to save the output JSON data opt.ArrayToStruct[0|1]: when set to 0, saveubjson outputs 1D/2D array in JSON array format; if sets to 1, an array will be shown as a struct with fields "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for sparse arrays, the non-zero elements will be saved to _ArrayData_ field in triplet-format i.e. (ix,iy,val) and "_ArrayIsSparse_" will be added with a value of 1; for a complex array, the _ArrayData_ array will include two columns (4 for sparse) to record the real and imaginary parts, and also "_ArrayIsComplex_":1 is added. opt.ParseLogical [1|0]: if this is set to 1, logical array elem will use true/false rather than 1/0. opt.SingletArray [0|1]: if this is set to 1, arrays with a single numerical element will be shown without a square bracket, unless it is the root object; if 0, square brackets are forced for any numerical arrays. opt.SingletCell [1|0]: if 1, always enclose a cell with "[]" even it has only one element; if 0, brackets are ignored when a cell has only 1 element. opt.ForceRootName [0|1]: when set to 1 and rootname is empty, saveubjson will use the name of the passed obj variable as the root object name; if obj is an expression and does not have a name, 'root' will be used; if this is set to 0 and rootname is empty, the root level will be merged down to the lower level. opt.JSONP [''|string]: to generate a JSONP output (JSON with padding), for example, if opt.JSON='foo', the JSON data is wrapped inside a function call as 'foo(...);' opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson back to the string form opt can be replaced by a list of ('param',value) pairs. The param string is equivallent to a field in opt and is case sensitive. output: examples: jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 'MeshTetra',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],... 'MeshTri',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;... 2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],... 'MeshCreator','FangQ','MeshTitle','T6 Cube',... 'SpecialData',[nan, inf, -inf]); saveubjson('jsonmesh',jsonmesh) saveubjson('jsonmesh',jsonmesh,'meshdata.ubj') license: BSD License, see LICENSE_BSD.txt files for details; implemented by `if(nargin==1)`.
- Lines 127-128: save to a file if FileName is set, suggested by Patrick Rapin; implemented by `filename=jsonopt('FileName','',opt)`.

### Key state/data transformations

- Lines 81: computes `varname` using `varname=inputname(1)`.
- Lines 82: computes `obj` using `obj=rootname`.
- Lines 86: computes `rootname` using `rootname=varname`.
- Lines 91: computes `opt` using `opt=struct('filename',varargin{1})`.
- Lines 95: computes `opt.IsOctave` using `opt.IsOctave=exist('OCTAVE_VERSION','builtin')`.
- Lines 99: computes `opt.singletarray` using `opt.singletarray=not(opt.norowbracket)`.
- Lines 102: computes `rootisarray` using `rootisarray=0`.
- Lines 103: computes `rootlevel` using `rootlevel=1`.
- Lines 104: computes `forceroot` using `forceroot=jsonopt('ForceRootName',0,opt)`.
- Lines 117: computes `json` using `json=obj2ubjson(rootname,obj,rootlevel,opt)`.
- Lines 122: computes `jsonp` using `jsonp=jsonopt('JSONP','',opt)`.
- Lines 128: computes `filename` using `filename=jsonopt('FileName','',opt)`.
- Lines 130: computes `fid` using `fid = fopen(filename, 'wb')`.

### Local helper functions

- Line 136: `obj2ubjson()` — `function txt=obj2ubjson(name,item,level,varargin)`.
  - Representative operation: `if(iscell(item))`.
  - Representative operation: `txt=cell2ubjson(name,item,level,varargin{:})`.
- Line 151: `cell2ubjson()` — `function txt=cell2ubjson(name,item,level,varargin)`.
  - Representative operation: `txt=''`.
  - Representative operation: `if(~iscell(item))`.
- Line 193: `struct2ubjson()` — `function txt=struct2ubjson(name,item,level,varargin)`.
  - Representative operation: `txt=''`.
  - Representative operation: `if(~isstruct(item))`.
- Line 243: `str2ubjson()` — `function txt=str2ubjson(name,item,level,varargin)`.
  - Representative operation: `txt=''`.
  - Representative operation: `if(~ischar(item))`.
- Line 277: `mat2ubjson()` — `function txt=mat2ubjson(name,item,level,varargin)`.
  - Representative operation: `if(~isnumeric(item) && ~islogical(item))`.
  - Representative operation: `error('input is not an array')`.
- Line 347: `matlabobject2ubjson()` — `function txt=matlabobject2ubjson(name,item,level,varargin)`. "st = struct(item);" would produce an inmutable warning, because it make the protected and private properties visible. Instead we get the
  - Representative operation: `if numel(item) == 0`.
  - Representative operation: `st = struct()`.
- Line 364: `matdata2ubjson()` — `function txt=matdata2ubjson(mat,level,varargin)`.
  - Representative operation: `if(isempty(mat))`.
  - Representative operation: `txt='Z'`.
- Line 415: `checkname()` — `function newname=checkname(name,varargin)`.
  - Representative operation: `isunpack=jsonopt('UnpackHex',1,varargin{:})`.
  - Representative operation: `newname=name`.

## Implementation structure

- json=saveubjson(rootname,obj,filename)
- json=saveubjson(rootname,obj,opt)
- json=saveubjson(rootname,obj,'param1',value1,'param2',value2,...)
- convert a MATLAB object (cell, struct or array) into a Universal
- Binary JSON (UBJSON) binary string
- author: Qianqian Fang (q.fang <at> neu.edu)
- created on 2013/08/17
- $Id$
- input:
- rootname: the name of the root-object, when set to '', the root name
- is ignored, however, when opt.ForceRootName is set to 1 (see below),
- the MATLAB variable name will be used as the root name.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `inputname()`, `ischar()`, `varargin2struct()`, `exist()`, `isfield()`, `not()`, `jsonopt()`, `islogical()`, `isstruct()`, `iscell()`, `isobject()`, `obj2ubjson()`, `fopen()`, `fwrite()`, `fclose()`, `cell2ubjson()`.
