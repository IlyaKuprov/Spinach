# interfaces/jsonlab-1.5/saveubjson.m

- Signature: `json=saveubjson(rootname,obj,varargin)`

## Purpose

No descriptive header was found. The best immediate identifier is `json=saveubjson(rootname,obj,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

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
