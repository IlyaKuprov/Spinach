# interfaces/jsonlab-1.5/mergestruct.m

- Signature: `s=mergestruct(s1,s2)`

## Purpose

Merge two scalar MATLAB structs, combining their fields.

## Physical / mathematical content

- The result starts with `s1`; fields from `s2` are then assigned, so `s2` values replace matching fields from `s1`.

## Numerical / algorithmic content

## Implementation structure

- s=mergestruct(s1,s2)
- merge two struct objects into one
- authors:Qianqian Fang (q.fang <at> neu.edu)
- date: 2012/12/22
- input:
- s1,s2: a struct object, s1 and s2 can not be arrays
- output:
- s: the merged struct object. fields in s1 and s2 will be combined in s.
- license:
- BSD License, see LICENSE_BSD.txt files for details
- --this function is part of jsonlab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
