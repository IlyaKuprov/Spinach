# kernel/eigenfields/rootmatch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/eigenfields/rootmatch.m`
- Signature: `[idx1,idx2,idx3]=rootmatch(field1,field2,field3,...`
- Total lines: 214

## Purpose

Global order-preserving root matching between three magnetic field root lists. The routine returns indices into the input lists that identify matching roots with minimum field-continuation cost. Syntax: [idx1,idx2,idx3]=rootmatch(field1,field2,field3,... edge12,edge23,edge31)

## Physical / mathematical content

- Eigenfield utilities. These files analyse field-dependent eigenstructure and resonance conditions, linking Hamiltonian spectra to magnetic-field sweeps and transition behaviour.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(field1,field2,field3,edge12,edge23,edge31)`.
- Lines 43-44: Initialise empty output; implemented by `idx1=[]; idx2=[]; idx3=[]`.
- Lines 46-47: Nothing to match; implemented by `if isempty(field1)||isempty(field2)||isempty(field3)`.
- Lines 51-52: Sort roots by field; implemented by `[field1,ord1]=sort(field1(:))`.
- Lines 56-57: Equal lists match by root order; implemented by `if (numel(field1)==numel(field2))&&(numel(field2)==numel(field3))`.
- Lines 64-65: Get list dimensions; implemented by `n1=numel(field1); n2=numel(field2); n3=numel(field3)`.
- Lines 67-68: Initialise dynamic programme; implemented by `match_count=zeros(n1+1,n2+1,n3+1)`.
- Lines 73-74: Dynamic programme over sorted root lists; implemented by `for a=1:(n1+1)`.
- Lines 78-79: Skip the initial state; implemented by `if (a==1)&&(b==1)&&(c==1), continue; end`.
- Lines 81-82: Initialise the local optimum; implemented by `best_count=-Inf`.
- Lines 86-87: Skip a root from the first list; implemented by `if (a>1)&&isfinite(match_cost(a-1,b,c))`.
- Lines 98-99: Skip a root from the second list; implemented by `if (b>1)&&isfinite(match_cost(a,b-1,c))`.
- Lines 110-111: Skip a root from the third list; implemented by `if (c>1)&&isfinite(match_cost(a,b,c-1))`.
- Lines 122-124: Match one root from each list; implemented by `if (a>1)&&(b>1)&&(c>1)&& isfinite(match_cost(a-1,b-1,c-1))`.
- Lines 138-139: Store the local optimum; implemented by `match_count(a,b,c)=best_count`.
- Lines 147-148: Start traceback; implemented by `a=n1+1; b=n2+1; c=n3+1`.
- Lines 150-151: Trace the optimal path; implemented by `while (a>1)||(b>1)||(c>1)`.
- Lines 153-154: Follow the stored move; implemented by `switch prev_move(a,b,c)`.

### Control flow inferred from the code

- Line 47: conditional branch on `isempty(field1)||isempty(field2)||isempty(field3)`.
- Line 57: conditional branch on `(numel(field1)==numel(field2))&&(numel(field2)==numel(field3))`.
- Line 74: `for` loop over `a=1:(n1+1)`.
- Line 75: `for` loop over `b=1:(n2+1)`.
- Line 76: `for` loop over `c=1:(n3+1)`.
- Line 79: conditional branch on `(a==1)&&(b==1)&&(c==1), continue; end`.
- Line 87: conditional branch on `(a>1)&&isfinite(match_cost(a-1,b,c))`.
- Line 90: conditional branch on `(test_count>best_count)||`.
- Line 99: conditional branch on `(b>1)&&isfinite(match_cost(a,b-1,c))`.
- Line 102: conditional branch on `(test_count>best_count)||`.
- Line 111: conditional branch on `(c>1)&&isfinite(match_cost(a,b,c-1))`.
- Line 114: conditional branch on `(test_count>best_count)||`.
- Line 123: conditional branch on `(a>1)&&(b>1)&&(c>1)&&`.
- Line 130: conditional branch on `(test_count>best_count)||`.

### Key state/data transformations

- Lines 44: computes `idx1` using `idx1=[]; idx2=[]; idx3=[]`.
- Lines 52: computes `[field1,ord1]` using `[field1,ord1]=sort(field1(:))`.
- Lines 53: computes `[field2,ord2]` using `[field2,ord2]=sort(field2(:))`.
- Lines 54: computes `[field3,ord3]` using `[field3,ord3]=sort(field3(:))`.
- Lines 59: computes `idx2` using `idx2=ord2.'`.
- Lines 60: computes `idx3` using `idx3=ord3.'`.
- Lines 65: computes `n1` using `n1=numel(field1); n2=numel(field2); n3=numel(field3)`.
- Lines 68: computes `match_count` using `match_count=zeros(n1+1,n2+1,n3+1)`.
- Lines 69: computes `match_cost` using `match_cost=inf(n1+1,n2+1,n3+1)`.
- Lines 70: computes `prev_move` using `prev_move=zeros(n1+1,n2+1,n3+1,'uint8')`.
- Lines 71: computes `match_cost(1,1,1)` using `match_cost(1,1,1)=0`.
- Lines 82: computes `best_count` using `best_count=-Inf`.
- Lines 83: computes `best_cost` using `best_cost=Inf`.
- Lines 84: computes `best_move` using `best_move=uint8(0)`.
- Lines 88: computes `test_count` using `test_count=match_count(a-1,b,c)`.
- Lines 89: computes `test_cost` using `test_cost=match_cost(a-1,b,c)`.
- Lines 125-127: computes `sheet_cost` using `sheet_cost=((field1(a-1)-field2(b-1))/edge12)^2+ ((field2(b-1)-field3(c-1))/edge23)^2+ ((field3(c-1)-field1(a-1))/edge31)^2`.
- Lines 139: computes `match_count(a,b,c)` using `match_count(a,b,c)=best_count`.

### Local helper functions

- Line 180: `grumble()` — `function grumble(field1,field2,field3,edge12,edge23,edge31)`.
  - Representative operation: `if (~isnumeric(field1))||(~isreal(field1))|| ((~isempty(field1))&&(~isvector(field1)))|| any(~isfinite(field1(:)))`.
  - Representative operation: `((~isempty(field1))&&(~isvector(field1)))|| any(~isfinite(field1(:)))`.

## Parameters / inputs

- field1 -real vector of roots at the first triangle vertex
- field2 -real vector of roots at the second triangle vertex
- field3 -real vector of roots at the third triangle vertex
- edge12 -positive distance between the first and the second
- triangle vertices
- edge23 -positive distance between the second and the third
- triangle vertices
- edge31 -positive distance between the third and the first
- triangle vertices

## Outputs

- idx1 -indices of matched roots in field1
- idx2 -indices of matched roots in field2
- idx3 -indices of matched roots in field3

## Implementation structure

- Global order-preserving root matching between three magnetic field root
- lists. The routine returns indices into the input lists that identify
- matching roots with minimum field-continuation cost. Syntax:
- [idx1,idx2,idx3]=rootmatch(field1,field2,field3,...
- edge12,edge23,edge31)
- field1 -real vector of roots at the first triangle vertex
- field2 -real vector of roots at the second triangle vertex
- field3 -real vector of roots at the third triangle vertex
- edge12 -positive distance between the first and the second
- triangle vertices
- edge23 -positive distance between the second and the third
- edge31 -positive distance between the third and the first

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `field1()`, `field2()`, `field3()`, `inf()`, `match_cost()`, `uint8()`, `match_count()`, `prev_move()`, `idx1()`, `ord1()`, `idx2()`, `ord2()`, `idx3()`, `ord3()`, `fliplr()`.
