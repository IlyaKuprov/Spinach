# interfaces/jsonlab-1.5/struct2jdata.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/jsonlab-1.5/struct2jdata.m`
- Signature: `newdata=struct2jdata(data,varargin)`
- Total lines: 96

## Purpose

No descriptive header was found. The best immediate identifier is `newdata=struct2jdata(data,varargin)`, and the implementation details below should be used to infer its role.

## Physical / mathematical content

- JSONLab vendored utilities. The main content is data serialisation, structure walking, option parsing, and text/binary JSON handling rather than spin physics.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 2-36: newdata=struct2jdata(data,opt,...) convert a JData object (in the form of a struct array) into an array input: data: a struct array. If data contains JData keywords in the first level children, these fields are parsed and regrouped into a data object (arrays, trees, graphs etc) based on JData specification. The JData keywords are "_ArrayType_", "_ArraySize_", "_ArrayData_" "_ArrayIsSparse_", "_ArrayIsComplex_" opt: (optional) a list of 'Param',value pairs for additional options The supported options include 'Recursive', if set to 1, will apply the conversion to every child; 0 to disable output: newdata: the covnerted data if the input data does contain a JData structure; otherwise, the same as the input. examples: obj=struct('_ArrayType_','double','_ArraySize_',[2 3], '_ArrayIsSparse_',1 ,'_ArrayData_',null); ubjdata=struct2jdata(obj); license: BSD License, see LICENSE_BSD.txt files for details; implemented by `fn=fieldnames(data)`.
- Lines 66-67: All-zeros sparse; implemented by `ndata=sparse(dim(1),prod(dim(2:end)))`.
- Lines 69-70: Sparse row vector; implemented by `ndata=sparse(1,ndata(:,1),ndata(:,2),dim(1),prod(dim(2:end)))`.
- Lines 72-73: Sparse column vector; implemented by `ndata=sparse(ndata(:,1),1,ndata(:,2),dim(1),prod(dim(2:end)))`.
- Lines 75-76: Generic sparse array.; implemented by `ndata=sparse(ndata(:,1),ndata(:,2),ndata(:,3),dim(1),prod(dim(2:end)))`.

### Control flow inferred from the code

- Line 40: `for` loop over `i=1:length(fn)`.
- Line 41: `for` loop over `j=1:len`.
- Line 50: `for` loop over `j=1:len`.
- Line 65: conditional branch on `isempty(ndata)`.

### Key state/data transformations

- Lines 36: computes `fn` using `fn=fieldnames(data)`.
- Lines 37: computes `newdata` using `newdata=data`.
- Lines 38: computes `len` using `len=length(data)`.
- Lines 43: computes `newdata(j)` using `newdata(j)=setfield(newdata(j),fn{i},jstruct2array(getfield(data(j),fn{i})))`.
- Lines 51: computes `ndata` using `ndata=cast(data(j).x0x5F_ArrayData_,data(j).x0x5F_ArrayType_)`.
- Lines 52: computes `iscpx` using `iscpx=0`.
- Lines 61: computes `dim` using `dim=double(data(j).x0x5F_ArraySize_)`.
- Lines 63: computes `ndata(:,end-1)` using `ndata(:,end-1)=complex(ndata(:,end-1),ndata(:,end))`.
- Lines 80: computes `ndata(:,3)` using `ndata(:,3)=complex(ndata(:,3),ndata(:,4))`.
- Lines 91: computes `newdata{j}` using `newdata{j}=ndata`.

## Implementation structure

- newdata=struct2jdata(data,opt,...)
- convert a JData object (in the form of a struct array) into an array
- authors:Qianqian Fang (q.fang <at> neu.edu)
- input:
- data: a struct array. If data contains JData keywords in the first
- level children, these fields are parsed and regrouped into a
- data object (arrays, trees, graphs etc) based on JData
- specification. The JData keywords are
- "_ArrayType_", "_ArraySize_", "_ArrayData_"
- "_ArrayIsSparse_", "_ArrayIsComplex_"
- opt: (optional) a list of 'Param',value pairs for additional options
- The supported options include

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fieldnames()`, `jsonopt()`, `isstruct()`, `getfield()`, `data()`, `newdata()`, `setfield()`, `jstruct2array()`, `strmatch()`, `cast()`, `double()`, `any()`, `ndata()`, `complex()`, `dim()`, `elseif()`.
