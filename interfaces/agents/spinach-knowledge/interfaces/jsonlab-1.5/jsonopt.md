# interfaces/jsonlab-1.5/jsonopt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/jsonopt.m`
- Signature: `val=jsonopt(key,default,varargin)`
- Total lines: 35

## Purpose

No descriptive header was found. The best immediate identifier is `val=jsonopt(key,default,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content

## Implementation structure

- val=jsonopt(key,default,optstruct)
- setting options based on a struct. The struct can be produced
- by varargin2struct from a list of 'param','value' pairs
- authors:Qianqian Fang (q.fang <at> neu.edu)
- $Id: loadjson.m 371 2012-06-20 12:43:06Z fangq $
- input:
- key: a string with which one look up a value from a struct
- default: if the key does not exist, return default
- optstruct: a struct where each sub-field is a key
- output:
- val: if key exists, val=optstruct.key; otherwise val=default
- license:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isstruct()`, `isfield()`, `getfield()`, `elseif()`, `lower()`.
