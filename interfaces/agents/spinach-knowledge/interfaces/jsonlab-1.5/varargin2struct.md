# interfaces/jsonlab-1.5/varargin2struct.m

- Signature: `opt=varargin2struct(varargin)`

## Purpose

Converts input arguments into one structure, combining name/value pairs with any struct arguments.

## Input and behavior

- `varargin` accepts parameter/value pairs and/or structures to merge. Each parameter name must be followed by its value.
- Parameter names supplied as strings are lowercased.
- With no inputs, the function returns an empty structure. Invalid input raises an error.

## Output

`opt` is the merged structure.

## Source notes

JSONLab utility by Qianqian Fang; created 2012/12/22. Distributed under the BSD License (see `LICENSE_BSD.txt`).
