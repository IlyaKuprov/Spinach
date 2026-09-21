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
