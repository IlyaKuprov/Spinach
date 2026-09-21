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
