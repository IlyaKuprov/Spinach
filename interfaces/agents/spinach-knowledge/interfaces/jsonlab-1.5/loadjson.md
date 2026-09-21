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
